"""Game hex grid (a longitude-wrapping equirectangular hex lattice).

This is the *game* grid, deliberately decoupled from the physics sphere grid.
Each hex centre maps to a lat/lon; resampling from the truth grid is nearest-cell
(area-aware resampling can replace this later).  A proper icosahedral Goldberg
hex grid can drop in without touching the compiler interface.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class HexGrid:
    width: int
    height: int

    @property
    def n(self) -> int:
        return self.width * self.height

    def latlon(self) -> tuple[np.ndarray, np.ndarray]:
        """Return (lat, lon) arrays of shape (height, width)."""
        # offset rows by half a cell to emulate hex stagger
        cols = np.arange(self.width)
        rows = np.arange(self.height)
        lon = -180.0 + (cols + 0.5) / self.width * 360.0
        lat = 90.0 - (rows + 0.5) / self.height * 180.0
        LAT, LON = np.meshgrid(lat, lon, indexing="ij")
        stagger = (np.arange(self.height) % 2) * (360.0 / self.width / 2.0)
        LON = LON + stagger[:, None]
        LON = ((LON + 180.0) % 360.0) - 180.0
        return LAT, LON

    def neighbors(self, r: int, c: int) -> list[tuple[int, int]]:
        """Axial-ish hex neighbours with longitude wrap."""
        even = (r % 2 == 0)
        if even:
            deltas = [(-1, -1), (-1, 0), (0, -1), (0, 1), (1, -1), (1, 0)]
        else:
            deltas = [(-1, 0), (-1, 1), (0, -1), (0, 1), (1, 0), (1, 1)]
        out = []
        for dr, dc in deltas:
            nr, nc = r + dr, (c + dc) % self.width
            if 0 <= nr < self.height:
                out.append((nr, nc))
        return out
