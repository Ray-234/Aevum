"""Drainage hydrology on the spherical graph.

Priority-Flood depression filling (so every land cell drains), steepest-descent
receivers, and flow accumulation of runoff.  Endorheic basins surface naturally
as cells whose filled level was raised above their bedrock.
"""
from __future__ import annotations

import heapq

import numpy as np


def priority_flood(grid, elev: np.ndarray, ocean: np.ndarray):
    """Return (filled, receiver, order_desc).

    ``filled`` is a hydrologically-correct surface; ``receiver[i]`` is the cell
    that ``i`` drains to (itself if it is an ocean sink); ``order_desc`` is a
    processing order from high to low filled elevation.
    """
    n = grid.n
    neigh = grid.neighbors
    filled = elev.copy()
    closed = np.zeros(n, dtype=bool)
    receiver = np.arange(n)
    heap: list[tuple[float, int]] = []

    for c in np.where(ocean)[0]:
        closed[c] = True
        heapq.heappush(heap, (filled[c], int(c)))

    if not heap:  # no ocean (e.g. fully frozen / dry): seed global minimum
        c = int(np.argmin(elev))
        closed[c] = True
        heapq.heappush(heap, (filled[c], c))

    order: list[int] = []
    while heap:
        level, c = heapq.heappop(heap)
        order.append(c)
        for nb in neigh[c]:
            if closed[nb]:
                continue
            closed[nb] = True
            if filled[nb] <= level:
                filled[nb] = level + 1e-3   # raise to spill level
            receiver[nb] = c
            heapq.heappush(heap, (filled[nb], int(nb)))
    order.reverse()                          # high -> low
    return filled, receiver, np.asarray(order, dtype=np.int64)


def flow_accumulation(grid, receiver: np.ndarray, order_desc: np.ndarray,
                      local_input: np.ndarray) -> np.ndarray:
    """Accumulate ``local_input`` (e.g. runoff volume) downstream."""
    acc = local_input.astype(np.float64).copy()
    for c in order_desc:
        r = receiver[c]
        if r != c:
            acc[r] += acc[c]
    return acc
