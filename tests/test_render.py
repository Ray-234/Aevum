import numpy as np
import pytest

from aevum import render


def test_elevation_legend_expands_lowland_bands():
    values = np.asarray([0.0, 200.0, 500.0, 1000.0, 2000.0, 3000.0, 4500.0])
    positions = np.asarray(render.ELEVATION_NORM(values), dtype=np.float64)

    assert render.ELEVATION_COLOR_TICKS == [
        -6000, -4500, -3000, -1500, 0, 200, 500, 1000, 2000, 3000, 4500, 6000
    ]
    assert np.all(np.diff(positions) > 0.06)
    assert positions[0] == pytest.approx(0.45)
    assert positions[1] == pytest.approx(0.53)
    assert positions[2] == pytest.approx(0.61)
