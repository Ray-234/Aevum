"""Spherical truth-layer grid (a pragmatic DGGS).

We use a spherical Fibonacci lattice for near-equal-area, hierarchy-friendly cell
centres, with neighbour adjacency from the convex hull (= spherical Delaunay) and
cell areas from a spherical Voronoi tessellation.  This is deliberately *not* the
game hex grid: physics runs here, the hex map is compiled from it via
area-conserving resampling, so the two can evolve resolution independently.

Three resolutions are intended (cf. project plan):
  * global evolution grid (~2e4 .. 8e4 cells)
  * local refinement grid (on demand)  -- not yet implemented
  * game hex grid (compiler.hexgrid)   -- separate module
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.spatial import ConvexHull, SphericalVoronoi, cKDTree


def fibonacci_sphere(n: int) -> np.ndarray:
    """Return ``(n, 3)`` near-uniform unit vectors on the sphere."""
    i = np.arange(n, dtype=np.float64) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)               # polar angle, equal-area in z
    golden = np.pi * (1.0 + 5.0 ** 0.5)              # golden angle
    theta = golden * i                               # azimuth
    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    return np.column_stack([x, y, z])


@dataclass
class SphereGrid:
    """An unstructured spherical grid with neighbour and area metadata."""

    xyz: np.ndarray          # (N, 3) unit vectors
    radius_m: float          # planet radius (sets physical areas / distances)

    # -- construction ----------------------------------------------------
    @classmethod
    def fibonacci(cls, n_cells: int, radius_m: float) -> "SphereGrid":
        return cls(xyz=fibonacci_sphere(n_cells), radius_m=float(radius_m))

    @property
    def n(self) -> int:
        return self.xyz.shape[0]

    # -- geographic coordinates -----------------------------------------
    @cached_property
    def lat(self) -> np.ndarray:
        return np.degrees(np.arcsin(np.clip(self.xyz[:, 2], -1.0, 1.0)))

    @cached_property
    def lon(self) -> np.ndarray:
        return np.degrees(np.arctan2(self.xyz[:, 1], self.xyz[:, 0]))

    @cached_property
    def latlon(self) -> np.ndarray:
        return np.column_stack([self.lat, self.lon])

    # -- adjacency (spherical Delaunay via convex hull) ------------------
    @cached_property
    def _hull(self) -> ConvexHull:
        return ConvexHull(self.xyz)

    @cached_property
    def neighbors(self) -> list[np.ndarray]:
        """List of neighbour index arrays, one per cell."""
        adj: list[set[int]] = [set() for _ in range(self.n)]
        for tri in self._hull.simplices:
            a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
            adj[a].update((b, c))
            adj[b].update((a, c))
            adj[c].update((a, b))
        return [np.array(sorted(s), dtype=np.int64) for s in adj]

    @cached_property
    def edges(self) -> np.ndarray:
        """Unique undirected edges as ``(E, 2)`` index pairs (i < j)."""
        s: set[tuple[int, int]] = set()
        for tri in self._hull.simplices:
            a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
            for u, v in ((a, b), (b, c), (a, c)):
                s.add((u, v) if u < v else (v, u))
        return np.array(sorted(s), dtype=np.int64)

    @cached_property
    def edge_lengths(self) -> np.ndarray:
        """Great-circle length (m) of each edge in :attr:`edges`."""
        i, j = self.edges[:, 0], self.edges[:, 1]
        dots = np.clip(np.einsum("ij,ij->i", self.xyz[i], self.xyz[j]), -1.0, 1.0)
        return np.arccos(dots) * self.radius_m

    # -- areas (spherical Voronoi) ---------------------------------------
    @cached_property
    def cell_area(self) -> np.ndarray:
        """Physical area (m^2) of each cell's Voronoi region."""
        sv = SphericalVoronoi(self.xyz, radius=1.0, center=np.zeros(3))
        solid_angle = sv.calculate_areas()           # steradians on unit sphere
        return solid_angle * (self.radius_m ** 2)

    @property
    def total_area(self) -> float:
        return float(self.cell_area.sum())

    # -- spatial query ---------------------------------------------------
    @cached_property
    def _kdtree(self) -> cKDTree:
        return cKDTree(self.xyz)

    def nearest(self, xyz: np.ndarray) -> np.ndarray:
        """Nearest cell index for each query unit vector ``(M, 3)``."""
        _, idx = self._kdtree.query(np.atleast_2d(xyz))
        return idx

    def nearest_latlon(self, lat: np.ndarray, lon: np.ndarray) -> np.ndarray:
        lat = np.radians(np.asarray(lat))
        lon = np.radians(np.asarray(lon))
        q = np.column_stack([
            np.cos(lat) * np.cos(lon),
            np.cos(lat) * np.sin(lon),
            np.sin(lat),
        ])
        return self.nearest(q)

    def great_circle_distance(self, i: int, j: int) -> float:
        d = float(np.clip(self.xyz[i] @ self.xyz[j], -1.0, 1.0))
        return np.arccos(d) * self.radius_m

    # -- differential operators (graph based) ---------------------------
    @cached_property
    def laplacian_weights(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return ``(edges, w)`` plus per-cell weight sum for an FV-style
        graph Laplacian.  Edge weight ~ 1/distance (cotangent-free, robust)."""
        w = 1.0 / np.maximum(self.edge_lengths, 1.0)
        wsum = np.zeros(self.n)
        np.add.at(wsum, self.edges[:, 0], w)
        np.add.at(wsum, self.edges[:, 1], w)
        return self.edges, w, wsum

    def diffuse(self, field: np.ndarray, kappa: float, dt: float,
                substeps: int = 1) -> np.ndarray:
        """Explicit graph diffusion of a scalar field (used by climate/erosion).

        ``kappa`` has units of m^2/s, ``dt`` in seconds.  Sub-stepping keeps the
        explicit scheme stable for large coefficients.
        """
        edges, w, wsum = self.laplacian_weights
        i, j = edges[:, 0], edges[:, 1]
        out = field.astype(np.float64).copy()
        sub_dt = dt / max(1, substeps)
        # Stability guard for the explicit Euler step.
        denom = kappa * sub_dt * np.maximum(wsum.max(), 1e-30)
        scale = min(1.0, 0.4 / denom) if denom > 0 else 1.0
        eff_dt = sub_dt * scale
        steps = max(1, int(np.ceil(substeps / scale))) if scale < 1 else substeps
        for _ in range(steps):
            flux = w * (out[j] - out[i])
            delta = np.zeros_like(out)
            np.add.at(delta, i, flux)
            np.add.at(delta, j, -flux)
            out += kappa * eff_dt * delta
        return out

    # -- resampling (area aware) ----------------------------------------
    def resample_to(self, other: "SphereGrid", field: np.ndarray) -> np.ndarray:
        """Nearest-cell resample of a scalar field onto ``other`` grid."""
        idx = self.nearest(other.xyz)
        return field[idx]
