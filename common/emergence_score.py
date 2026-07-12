"""
common/emergence_score.py — validated fast screening score, generalized 2D/3D.

PURPOSE: a fast, discriminating score of HOW FAR a run self-organizes from t=0, to RANK candidate
initial conditions during a mass-random search. Feeds the search ledger alongside
common/emergence.compute_level_report's reached_level (the two-stage design this repo's discipline
requires: this score is a cheap screen, reached_level is the audited Level 0-8 measurement).

PROVENANCE: ported from a 2D-validated score (emergence_score.py, supplied by Claude's sandbox
work) with the exact same logic (entropy/sf_peak/blobs/churn/L1-L4/complexity-window/sweet-spot/
penalties). The only change here is dimension-agnostics: `_coarse_entropy` block-averages over
however many spatial axes the field has (2D: ny,nx block=2D; 3D: nz,ny,nx block=3D cubes),
`_sf_peak` computes the FFT power spectrum and radial |k| binning over ND axes (uses
common.diagnostics._radial_power_spectrum, already ND per its docstring), and the blob count uses
scipy.ndimage.label which is natively ND. The 2D self-test in this module's __main__ reproduces the
original validated numbers so this is a regression-checked generalization, not a rewrite.

HONEST FLOOR: this is a fast SCREEN / proxy that ranks candidates -- NOT the full Emergence Level
0-8 with tracked individuality. Deep Levels (persistent-individuality with identity tracking,
self-maintaining closure, division/inheritance) require Room-specific object tracking
(common/emergence.py + per-Room trackers). Use this score to RANK which ICs deserve deeper
evaluation; keep BOTH this score and reached_level in the ledger.

USAGE:
    from common.emergence_score import emergence_score
    traj = [snap_t0, snap_mid, snap_final]   # >=2 fields (2D or 3D, real or complex), time-ascending
    total, detail = emergence_score(traj)
"""

import os
import sys

import numpy as np
from scipy.ndimage import label

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common.diagnostics import _radial_power_spectrum  # noqa: E402


def _coarse_entropy(amp, block=4):
    """Coarse-grained spatial entropy in [0,1], ND-agnostic (2D: block x block; 3D: block^3 cubes).
    Uniform->0 (guarded), structured->~0.6-0.8, noise->~1. Identical logic/thresholds to the
    validated 2D version, just block-averaging over however many axes `amp` has."""
    shape = amp.shape
    m = [n // block for n in shape]
    if any(mi < 1 for mi in m):
        return 0.0
    trimmed = amp[tuple(slice(0, mi * block) for mi in m)]
    # reshape (m0,block,m1,block,...) then average over all the "block" axes (odd positions)
    reshaped = trimmed.reshape([v for mi in m for v in (mi, block)])
    block_axes = tuple(range(1, 2 * len(m), 2))
    coarse = reshaped.mean(axis=block_axes)
    contrast = (coarse.max() - coarse.min()) / (coarse.mean() + 1e-9)
    if contrast < 0.05:
        return 0.0
    c = coarse - coarse.min()
    p = c / c.sum()
    p = p[p > 0]
    return float(-(p * np.log(p)).sum() / np.log(len(p)))


def _sf_peak(amp):
    """Structure-factor peak prominence (ND): high when there is a dominant length scale.
    Reuses common.diagnostics._radial_power_spectrum (already ND, shell-binned by integer |k|)."""
    kbins, radial, counts = _radial_power_spectrum(amp)
    sf = radial / np.maximum(counts, 1)
    if len(sf) < 20:
        return 0.0
    return float(sf[1:15].max() / (sf[1:25].mean() + 1e-9))


def emergence_score(traj,
                    weights=(1.0, 1.2, 1.5, 2.0, 1.0),   # (L1, L2, L3, L4, window)
                    sweet_spot=(2, 12), sweet_bonus=1.0,
                    churn_ceiling=0.25, turbulent_penalty=0.5,
                    frag_ceiling=40, frag_penalty=0.5):
    """
    Score how far a run self-organized. traj = list of >=2 ND (2D or 3D) snapshots (real or complex;
    last two used for churn). Returns (total_score: float, detail: dict). Higher = closer to
    'seed -> plant' (persistent structured individuals in the complexity window). Dead/exploded/
    turbulent/uniform -> low. Logic identical to the validated 2D score; only the entropy/sf_peak/
    blob primitives are ND-generalized (see module docstring).
    """
    psiF = traj[-1]
    amp = np.abs(psiF)
    mean = amp.mean()
    d = {}
    if not (np.isfinite(mean) and 0.02 < mean < 50):
        d["status"] = "dead" if (np.isfinite(mean) and mean <= 0.02) else "exploded/NaN"
        return 0.0, d
    rel_var = float(amp.var() / (mean ** 2 + 1e-9))
    prom = _sf_peak(amp)
    churn = float(np.abs(np.abs(traj[-1]) - np.abs(traj[-2])).mean() / (mean + 1e-9))
    blobs = int(label(amp > 0.5 * amp.max())[1])
    ent = _coarse_entropy(amp)

    L1 = rel_var > 0.02
    L2 = prom > 2.5
    L3 = (churn < churn_ceiling) and L2
    L4 = (2 <= blobs <= 30) and L2 and L3
    in_window = (0.5 < ent < 0.97) and L2

    w = weights
    s = w[0] * L1 + w[1] * L2 + w[2] * L3 + w[3] * L4 + w[4] * in_window
    if sweet_spot[0] <= blobs <= sweet_spot[1] and L3 and in_window:
        s += sweet_bonus
        d["sweet_spot"] = True
    if churn > churn_ceiling:
        s -= turbulent_penalty
        d["turbulent"] = True
    if blobs > frag_ceiling:
        s -= frag_penalty
        d["fragmented"] = True

    d.update(rel_var=round(rel_var, 3), sf_peak=round(prom, 2), churn=round(churn, 3),
             blobs=blobs, coarse_entropy=round(ent, 3), window=bool(in_window),
             levels=dict(L1=bool(L1), L2=bool(L2), L3=bool(L3), L4=bool(L4)))
    return round(max(s, 0.0), 2), d


if __name__ == "__main__":
    # Regression check against the original 2D-validated numbers, then a 3D analog.
    print("--- 2D regression check (should match the validated 2D module) ---")
    N = 64
    X, Y = np.meshgrid(np.arange(N), np.arange(N), indexing="ij")
    rng = np.random.default_rng(0)
    few = sum(np.exp(-((X - cx) ** 2 + (Y - cy) ** 2) / 25) for cx, cy in [(16, 16), (16, 48), (48, 32)])
    few_traj = [(few + 0.01 * rng.standard_normal((N, N))).astype(complex) for _ in range(3)]
    dead_traj = [np.full((N, N), 0.001, complex)] * 3
    noise_traj = [(rng.standard_normal((N, N)) + 1j * rng.standard_normal((N, N))) for _ in range(3)]
    print("few_blobs (should be high ~7):", emergence_score(few_traj)[0])
    print("dead      (should be 0):      ", emergence_score(dead_traj)[0])
    print("noise     (should be low):    ", emergence_score(noise_traj)[0])

    print("\n--- 3D analog ---")
    N3 = 32
    X3, Y3, Z3 = np.meshgrid(np.arange(N3), np.arange(N3), np.arange(N3), indexing="ij")
    centers3 = [(8, 8, 8), (8, 24, 24), (24, 16, 8)]
    few3 = sum(np.exp(-((X3 - cx) ** 2 + (Y3 - cy) ** 2 + (Z3 - cz) ** 2) / 20) for cx, cy, cz in centers3)
    few3_traj = [(few3 + 0.01 * rng.standard_normal((N3, N3, N3))).astype(complex) for _ in range(3)]
    dead3_traj = [np.full((N3, N3, N3), 0.001, complex)] * 3
    noise3_traj = [(rng.standard_normal((N3, N3, N3)) + 1j * rng.standard_normal((N3, N3, N3))) for _ in range(3)]
    print("3D few_blobs (should be high):", emergence_score(few3_traj)[0])
    print("3D dead      (should be 0):   ", emergence_score(dead3_traj)[0])
    print("3D noise     (should be low): ", emergence_score(noise3_traj)[0])
