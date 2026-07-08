"""Process modules.  Each implements the unified ``Module.step`` contract.

First-version modules are intentionally low fidelity but causally complete:
every mountain, river, desert and ore body traces to a process and an event,
not to raw noise.  Modules can later be swapped for higher-fidelity kernels
without changing the architecture.
"""

from aevum.modules.stellar import StellarModule
from aevum.modules.interior import InteriorModule
from aevum.modules.tectonics import TectonicsModule
from aevum.modules.impacts import ImpactModule
from aevum.modules.terrain import TerrainModule
from aevum.modules.climate import ClimateModule
from aevum.modules.biogeochem import BiogeochemModule
from aevum.modules.biosphere import BiosphereModule
from aevum.modules.resources import ResourceModule

__all__ = [
    "StellarModule",
    "InteriorModule",
    "TectonicsModule",
    "ImpactModule",
    "TerrainModule",
    "ClimateModule",
    "BiogeochemModule",
    "BiosphereModule",
    "ResourceModule",
]
