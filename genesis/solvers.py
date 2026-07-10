"""Shared pseudo-spectral numerics for the G001 / G003 / dividing-protocell 3D Room solvers written
against `requests/Codex計算依頼書_Genesis3D探索_G003_G001_protocell.md`. These are solver-specific
helpers (Codex/Claude's role under docs/AGENTS.md, not the fixed `common/` scaffold).
"""

import numpy as np


def k_grid(shape, spacing=1.0):
    """Angular wavenumber component grids (2*pi*fftfreq/spacing, meshgrid 'ij') and k^2 = sum(k_i^2)."""
    ks = [2 * np.pi * np.fft.fftfreq(n, d=spacing) for n in shape]
    kk = np.meshgrid(*ks, indexing="ij")
    k2 = sum(k ** 2 for k in kk)
    return kk, k2


def dealias_mask(shape):
    """2/3-rule dealiasing mask: keep the inner 2/3 of each axis's wavenumber range (cycles/box)."""
    masks = []
    for n in shape:
        freqs = np.fft.fftfreq(n) * n
        masks.append(np.abs(freqs) <= n / 3.0)
    grids = np.meshgrid(*masks, indexing="ij")
    out = grids[0]
    for g in grids[1:]:
        out = out & g
    return out


def leray_project(u_hat, kk, k2):
    """Project a vector field (list of component FFTs) onto its divergence-free part (spectral
    Helmholtz-Leray projection): u_hat -= k*(k.u_hat)/k^2, with the k=0 (mean) mode left untouched.
    """
    k_dot_u = sum(kk[a] * u_hat[a] for a in range(len(kk)))
    safe_k2 = np.where(k2 == 0, 1.0, k2)
    factor = np.where(k2 == 0, 0.0, k_dot_u / safe_k2)
    return [u_hat[a] - kk[a] * factor for a in range(len(kk))]
