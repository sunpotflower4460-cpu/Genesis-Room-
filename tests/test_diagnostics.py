"""Minimal tests for common/diagnostics.py: verify each measurement against a hand-built field with
a known answer (not against a physics run -- these are unit tests of the measurement tools themselves).
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import diagnostics as diag  # noqa: E402


def test_variance_growth_saturating_complex_field():
    """Regression test: amplitude that grows away from a near-zero noise floor then SATURATES (e.g.
    CGL/GPE |A|->1) must report a positive growth_rate from the early departure, not a naive
    full-trajectory slope that goes negative once late-time coarsening/relaxation dominates the fit.
    """
    t = np.linspace(0, 10, 40)
    amp = np.tanh(t)  # grows 0->1 then flat; a real run's variance also DECLINES a bit past the peak
    decline = np.concatenate([amp[:25], amp[25] - 0.05 * (t[25:] - t[25])])
    rng = np.random.default_rng(0)
    fields = [decline[i] * np.exp(1j * rng.uniform(0, 2 * np.pi, (16, 16))) for i in range(len(t))]
    variances, growth_rate = diag.variance_growth(fields)
    assert growth_rate > 0
    assert variances[0] < variances[int(np.argmax(variances))]


def test_variance_growth_real_field_preserves_sign_structure():
    """A real +/-1 double-well field (e.g. Cahn-Hilliard phi) must use Var[field] directly -- Var[|field|]
    would erase the +/- domain structure that IS the Level-1 signal.
    """
    n = 32
    noise = 0.01 * np.ones((n, n))
    saturated = np.where(np.indices((n, n))[0] < n // 2, 1.0, -1.0)
    variances, growth_rate = diag.variance_growth([noise, saturated])
    assert variances[-1] > variances[0]
    assert growth_rate > 0


def test_structure_factor_peak_known_wavelength():
    n = 64
    wavelength = 8.0
    x = np.arange(n)
    field = np.tile(np.sin(2 * np.pi * x / wavelength)[:, None], (1, n))
    peak_k, prominence = diag.structure_factor_peak(field)
    expected_k = n / wavelength
    assert abs(peak_k - expected_k) <= 1.0
    assert prominence > 1.5


def _single_vortex_2d(n=48, sigma=6.0):
    center = (n / 2.0) - 0.5  # sits between grid points, so the singularity isn't ON a node
    y, x = np.meshgrid(np.arange(n), np.arange(n), indexing="ij")
    dy, dx = y - center, x - center
    r2 = dx ** 2 + dy ** 2
    amp = np.exp(-r2 / (2 * sigma ** 2))
    theta = np.arctan2(dy, dx)
    return (amp * np.exp(1j * theta)).astype(np.complex128)


def test_winding_defect_count_single_vortex_2d():
    field = _single_vortex_2d()
    assert diag.winding_defect_count(field) == 1


def test_winding_defect_count_straight_line_vortex_3d():
    nz = 5
    vortex_2d = _single_vortex_2d()
    field_3d = np.repeat(vortex_2d[None, :, :], nz, axis=0)
    assert diag.winding_defect_count(field_3d) == nz


def test_connected_components_known_count():
    mask = np.zeros((30, 30), dtype=bool)
    mask[2:5, 2:5] = True
    mask[10:14, 10:13] = True
    mask[20:25, 20:24] = True
    count, _, sizes = diag.connected_components(mask)
    assert count == 3
    assert sorted(sizes.tolist()) == sorted([9, 12, 20])


def test_coarsening_length_monotonic():
    rng = np.random.default_rng(0)
    n = 64

    fine = rng.choice([-1.0, 1.0], size=(n, n))

    block = 8
    coarse_blocks = rng.choice([-1.0, 1.0], size=(n // block, n // block))
    coarse = np.repeat(np.repeat(coarse_blocks, block, axis=0), block, axis=1)

    l_fine = diag.coarsening_length(fine)
    l_coarse = diag.coarsening_length(coarse)
    assert l_coarse > l_fine


if __name__ == "__main__":
    test_structure_factor_peak_known_wavelength()
    test_winding_defect_count_single_vortex_2d()
    test_winding_defect_count_straight_line_vortex_3d()
    test_connected_components_known_count()
    test_coarsening_length_monotonic()
    print("all tests passed")
