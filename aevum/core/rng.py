"""Deterministic, hierarchical random-number keying.

The random seed for any draw is derived from
``planet_seed + module_name + time + event_index``.  Because each module hashes
its own namespace, adding a *new* module later does not shift the random stream
of existing modules -- old worlds reproduce bit-for-bit except where the new
module actually intervenes.
"""
from __future__ import annotations

import hashlib
import struct
from dataclasses import dataclass

import numpy as np


def derive_seed(planet_seed: int, *parts: object) -> int:
    """Derive a stable 64-bit seed from a planet seed and arbitrary key parts."""
    h = hashlib.blake2b(digest_size=8)
    h.update(struct.pack("<q", int(planet_seed) & 0x7FFF_FFFF_FFFF_FFFF))
    for p in parts:
        if isinstance(p, float):
            # Quantise floats so that 1.0 and 1.0000000001 do not diverge.
            p = round(p, 6)
        h.update(b"\x00")
        h.update(str(p).encode("utf-8"))
    return int.from_bytes(h.digest(), "little")


def generator(planet_seed: int, module: str, time_myr: float = 0.0,
              event_index: int = 0) -> np.random.Generator:
    """Return a numpy Generator for a (module, time, event) namespace."""
    return np.random.default_rng(derive_seed(planet_seed, module, time_myr, event_index))


@dataclass(frozen=True)
class RNGKey:
    """An opaque, composable key handed to modules each step.

    A module derives sub-streams via :meth:`child` so that different sub-processes
    (e.g. rifting vs. subduction within tectonics) stay independent and stable.
    """

    planet_seed: int
    module: str
    time_myr: float
    event_index: int = 0

    def generator(self) -> np.random.Generator:
        return generator(self.planet_seed, self.module, self.time_myr, self.event_index)

    def child(self, label: str, index: int = 0) -> "RNGKey":
        return RNGKey(
            planet_seed=self.planet_seed,
            module=f"{self.module}/{label}",
            time_myr=self.time_myr,
            event_index=index,
        )
