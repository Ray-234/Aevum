"""Resource genesis with a full "formation -> preservation" lineage.

Deposits are objects, not ``copper = 1`` on a cell.  Each records its genesis
model, formation age, host lithology, geometry, commodity vector, grade/tonnage,
burial depth and preservation/erosion history.  "Modern resources" are therefore
a *query*:

    available = original_event x preservation x burial x exposure x technology

The same physical deposit can read as absent / known / exploitable / low-grade
depending on the civilisation's technology stage (handled by the compiler).
"""
from __future__ import annotations

import numpy as np

from aevum.core.events import Event
from aevum.core.module import Module, StepResult
from aevum.core.provenance import Provenance


class ResourceModule(Module):
    name = "resources"
    produces = ["resources.deposits"]
    fidelity = "genesis_rules"
    interval_myr = 50.0

    def init_state(self, world, rng_key) -> None:
        world.objects["resources.deposits"] = []
        self._dep_id = 0

    def _add(self, world, deposits, t, cell, model, commodities, grade, tonnage,
             host, events, cause):
        dep = {
            "id": self._dep_id,
            "genesis_model": model,
            "formation_age_myr": round(t, 1),
            "cell": int(cell),
            "host_lithology": float(host),
            "commodity_vector": commodities,
            "grade": round(float(grade), 3),
            "tonnage_index": round(float(tonnage), 2),
            "burial_depth_m": 0.0,
            "preservation": 1.0,
            "formation_sediment_m": float(world.get_field("sediment.thickness_m")[cell]),
            "formation_erosion_m": float(world.get_field("erosion_m", 0.0)[cell]),
        }
        self._dep_id += 1
        deposits.append(dep)
        events.append(Event("mineralization", t, self.name, location=int(cell),
                            magnitude=grade, params={"model": model,
                                                     "commodities": list(commodities),
                                                     "cause": cause}))
        return dep

    def step(self, world, t, dt, forcing, rng_key) -> StepResult:
        grid = world.grid
        rng = rng_key.generator()
        deposits = world.object_set("resources.deposits")

        litho = world.get_field("terrain.lithology", 0.0)
        volc_age = world.get_field("tectonics.volcanism_age_myr", -1.0)
        orog_age = world.get_field("tectonics.orogeny_age_myr", -1.0)
        sediment = world.get_field("sediment.thickness_m", 0.0)
        erosion = world.get_field("erosion_m", 0.0)
        accum = world.get_field("hydrology.flow_accumulation", 0.0)
        precip = world.get_field("climate.precipitation", 500.0)
        evap = world.get_field("climate.evaporation", 0.0)
        npp = world.get_field("biosphere.npp", 0.0)
        ocean = world.ocean_mask()
        land = ~ocean
        o2 = world.g("biogeochem.oxygen_fraction", 0.0)

        events: list[Event] = []
        active_arc = volc_age >= (t - 60.0)
        active_orog = orog_age >= (t - 80.0)

        def sample(mask, k):
            cells = np.where(mask)[0]
            if cells.size == 0:
                return []
            return rng.choice(cells, size=min(k, cells.size), replace=False)

        # hydrothermal porphyry/VMS at active arcs
        for c in sample(active_arc & land, 4):
            self._add(world, deposits, t, c, "hydrothermal_porphyry",
                      {"copper": 1.0, "gold": 0.3}, rng.uniform(0.3, 1.2),
                      rng.uniform(0.5, 3.0), litho[c], events,
                      "magmatic fluids above subduction arc")
        # magmatic Ni/Cr/PGE at igneous provinces
        for c in sample(active_arc & (litho == 0.0), 2):
            self._add(world, deposits, t, c, "magmatic_intrusion",
                      {"nickel": 1.0, "chromium": 0.6, "platinum": 0.2},
                      rng.uniform(0.2, 0.8), rng.uniform(0.4, 2.0), litho[c], events,
                      "fractional crystallisation in mafic intrusion")
        # banded iron: anoxic ocean once photosynthesis exists (pre/around GOE)
        if 0.0 < o2 < 0.02 and world.g("bio.oxygenic_photosynthesis", 0.0) > 0.5:
            for c in sample(ocean, 3):
                self._add(world, deposits, t, c, "banded_iron_formation",
                          {"iron": 1.0}, rng.uniform(0.4, 0.9), rng.uniform(2.0, 6.0),
                          1.0, events, "oxidation of dissolved Fe in low-O2 ocean")
        # evaporites in arid endorheic settings
        arid = land & (precip < 250) & (evap > precip + 200)
        for c in sample(arid, 3):
            self._add(world, deposits, t, c, "evaporite",
                      {"potash": 1.0, "salt": 0.8, "lithium": 0.2},
                      rng.uniform(0.3, 0.7), rng.uniform(0.5, 2.5), 1.0, events,
                      "evaporation of restricted basin")
        # placer gold downstream of orogens / igneous
        thresh = np.percentile(accum[land], 95) if land.any() else np.inf
        for c in sample(land & (accum > thresh) & (active_orog | (litho == 0.0)), 2):
            self._add(world, deposits, t, c, "placer",
                      {"gold": 0.7, "titanium": 0.3}, rng.uniform(0.1, 0.4),
                      rng.uniform(0.2, 1.0), 1.0, events,
                      "fluvial concentration of dense minerals")
        # coal: vascular forests on land with burial
        if world.g("bio.vascular_forests", 0.0) > 0.5:
            for c in sample(land & (npp > 0.6) & (sediment > 20), 3):
                self._add(world, deposits, t, c, "coal",
                          {"coal": 1.0}, rng.uniform(0.5, 1.0), rng.uniform(1.0, 4.0),
                          1.0, events, "burial of terrestrial plant biomass")
        # oil/gas: marine productivity + anoxia + burial
        if o2 > 0.05:
            src = ocean & (npp > 0.4) & (sediment > 30)
            for c in sample(src, 3):
                self._add(world, deposits, t, c, "petroleum_system",
                          {"oil": 1.0, "gas": 0.5}, rng.uniform(0.3, 0.9),
                          rng.uniform(1.0, 5.0), 1.0, events,
                          "burial+maturation of marine organic matter")

        # --- update preservation/burial of existing deposits ------------
        for dep in deposits:
            c = dep["cell"]
            buried = max(sediment[c] - dep["formation_sediment_m"], 0.0)
            exhumed = max(erosion[c] - dep["formation_erosion_m"], 0.0)
            dep["burial_depth_m"] = round(float(buried), 1)
            # erosion can destroy shallow/surface deposits
            if exhumed > 1500.0 and dep["genesis_model"] in ("placer", "evaporite"):
                dep["preservation"] = round(max(dep["preservation"] - 0.1, 0.0), 3)

        deposits[:] = [d for d in deposits if d["preservation"] > 0.05]

        world.provenance.record(Provenance(
            "resources.deposits", self.name, self.fidelity, "1",
            direct_cause=f"{len(events)} new deposits this step; {len(deposits)} extant",
            upstream_events=[e.id for e in events[:5]]))
        diag = {"n_new": len(events), "n_total": len(deposits)}
        return StepResult(state_delta={"objects": {"resources.deposits": deposits}},
                          events=events, diagnostics=diag)
