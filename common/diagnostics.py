"""Shared 2D/3D diagnostics for measuring docs/EMERGENCE_LEVELS.md indicators BY NUMBER, not by
looking at an image. Every function is a pure function of numpy array input -> numeric output(s).
No file I/O, no plotting, no side effects.

Field arrays are 2D (ny, nx) or 3D (nz, ny, nx) unless noted. Complex fields carry an order-parameter
phase/amplitude (e.g. psi for TDGL / CGL / GPE); real fields carry a scalar order parameter
(e.g. phi for Cahn-Hilliard / Model H).
"""

import numpy as np
from scipy import ndimage


def _wrap_angle(x):
    """Wrap an angle (difference) into (-pi, pi]."""
    return (x + np.pi) % (2 * np.pi) - np.pi


def _radial_power_spectrum(field):
    """FFT power spectrum of `field` (mean removed) binned onto integer radial wavenumber shells.
    Returns (kbin_values, radial_power, radial_counts) with kbin_values[0] == 0 (the DC shell).
    """
    f = np.asarray(field)
    spec = np.fft.fftn(f - f.mean())
    power = np.abs(spec) ** 2
    freqs = [np.fft.fftfreq(n) * n for n in f.shape]
    grids = np.meshgrid(*freqs, indexing="ij")
    kr = np.sqrt(sum(g ** 2 for g in grids))
    kbin = kr.round().astype(int)
    nbins = kbin.max() + 1
    radial = np.bincount(kbin.ravel(), weights=power.ravel(), minlength=nbins)
    counts = np.bincount(kbin.ravel(), minlength=nbins)
    return np.arange(nbins), radial, counts


def variance_growth(field_series):
    """Lv1: growth rate of spatial variance over time (EMERGENCE_LEVELS.md Level 1: "空間分散が
    ノイズ床を超えて増大" -- var_growth > 0 is part of the Level-1 gate `var_growth > 0 AND ...`).

    field_series: sequence of T field snapshots (each 2D or 3D, real or complex), earliest first.
    Returns (variances, growth_rate): per-frame spatial variance of |field| and the least-squares
    slope of variance vs frame index (positive = growing away from the noise floor).
    """
    variances = np.array([float(np.var(np.abs(np.asarray(f)))) for f in field_series])
    t = np.arange(len(variances))
    if len(t) < 2:
        return variances, 0.0
    growth_rate = float(np.polyfit(t, variances, 1)[0])
    return variances, growth_rate


def structure_factor_peak(field):
    """Lv1: structure factor S(k) peak wavenumber k* (EMERGENCE_LEVELS.md Level 1: "構造因子 S(k) に
    ピーク（特徴波長 k*）が立つ"). 2D/3D via FFT, radially binned by integer wavenumber shell.

    Returns (peak_k, prominence): the shell index (k=0 dropped) with maximum mean power, and its
    prominence = peak_power / mean_power_over_shells (>~1.5 indicates a real peak, not flat noise).
    """
    _, radial, counts = _radial_power_spectrum(field)
    nbins = radial.shape[0]
    if nbins <= 2:
        return 0.0, 0.0
    radial_mean = radial / np.maximum(counts, 1)
    profile = radial_mean[1:]  # drop k=0 (the mean)
    peak_idx = int(np.argmax(profile)) + 1
    mean_val = float(profile.mean()) if profile.size else 0.0
    prominence = float(radial_mean[peak_idx] / mean_val) if mean_val > 0 else 0.0
    return float(peak_idx), prominence


def correlation_length(field):
    """Lv1: correlation length xi (EMERGENCE_LEVELS.md Level 1: "相関長 ξ が有限/成長").
    xi = 1 / k*, the inverse of the structure-factor peak wavenumber. Returns 0.0 if S(k) is flat
    (no peak -> Level 0, undefined correlation length).
    """
    peak_k, _ = structure_factor_peak(field)
    return float(1.0 / peak_k) if peak_k > 0 else 0.0


def coarsening_length(field):
    """Lv2/Lv5: characteristic domain size L = 2*pi / <k> using the first moment of S(k)
    (EMERGENCE_LEVELS.md Level 2 "構造寿命" / Level 5 co-differentiated domains coarsening in time).
    Larger domains (coarser structure) -> smaller mean k -> larger L.
    """
    kvals, radial, counts = _radial_power_spectrum(field)
    freqs_flat = kvals  # shell index == radial |k| in cycles/box for this binning
    power = radial
    mask = freqs_flat > 0
    total_power = float(power[mask].sum())
    if not np.any(mask) or total_power <= 0:
        return 0.0
    mean_k = float(np.sum(freqs_flat[mask] * power[mask]) / total_power)
    return float(2 * np.pi / mean_k) if mean_k > 0 else 0.0


def _plaquette_winding(theta, axis_a, axis_b, amp=None, thr=0.0):
    """Integer plaquette winding number of phase array `theta` around the unit loop in the
    (axis_a, axis_b) plane, batched over all remaining axes. Returns an int array shaped like theta;
    winding[...] != 0 marks a plaquette pierced by a phase singularity (a point vortex in 2D, or one
    lattice cell of a vortex line's path in 3D).

    If `amp` (field amplitude) and `thr` are given, plaquettes whose mean corner amplitude is at or
    below `thr` are zeroed out -- this masks the disordered/noise region so it cannot register as a
    spurious defect (a vortex core is low-amplitude but its SURROUNDING plaquette sits in the ordered
    bulk; pure noise is not).
    """
    a = theta
    b = np.roll(theta, -1, axis_a)
    c = np.roll(np.roll(theta, -1, axis_a), -1, axis_b)
    d = np.roll(theta, -1, axis_b)
    circ = _wrap_angle(b - a) + _wrap_angle(c - b) + _wrap_angle(d - c) + _wrap_angle(a - d)
    winding = np.round(circ / (2 * np.pi)).astype(int)
    if amp is not None:
        mean4 = 0.25 * (amp + np.roll(amp, -1, axis_a)
                         + np.roll(np.roll(amp, -1, axis_a), -1, axis_b)
                         + np.roll(amp, -1, axis_b))
        winding = np.where(mean4 > thr, winding, 0)
    return winding


def winding_defect_count(complex_field, amplitude_threshold_frac=0.25):
    """Lv2: topological phase-winding defects, NOT mere low density (EMERGENCE_LEVELS.md Level 2:
    "欠陥数（巻き数・位相特異点。単なる低密度でなく位相巻きで検出）"). A plaquette whose net phase
    circulation is a nonzero multiple of 2*pi is pierced by a phase singularity.

    2D field (ny, nx): counts plaquettes with nonzero winding on the single (0,1) lattice -- point
    vortices.
    3D field (nz, ny, nx): a vortex LINE pierces one face per lattice cell along its path. Per this
    module's docstring contract, we compute plaquette winding on all three face orientations (xy:
    axes (1,2); yz: axes (0,1); zx: axes (0,2)) and sum the nonzero-winding faces -- a lower-bound
    proxy for total vortex-line length in lattice units (a line running straight along z pierces one
    (y,x)-face per z-layer it passes through). This is an assumption documented here because a single
    scalar "line count" requires an explicit line-tracing algorithm this module does not implement;
    the face-piercing count is verified against a known straight-line vortex in tests/test_diagnostics.py.

    Only plaquettes in the ordered bulk (mean corner amplitude above `amplitude_threshold_frac` of the
    field's peak amplitude) count, so disordered near-zero noise doesn't register as spurious defects.
    """
    field = np.asarray(complex_field)
    theta = np.angle(field)
    amp = np.abs(field)
    thr = amplitude_threshold_frac * float(amp.max()) if amp.size else 0.0
    if field.ndim == 2:
        w = _plaquette_winding(theta, 0, 1, amp, thr)
        return int(np.count_nonzero(w))
    if field.ndim == 3:
        total = 0
        for axis_a, axis_b in [(1, 2), (0, 1), (0, 2)]:
            w = _plaquette_winding(theta, axis_a, axis_b, amp, thr)
            total += int(np.count_nonzero(w))
        return total
    raise ValueError("winding_defect_count expects a 2D or 3D field, got ndim=%d" % field.ndim)


def connected_components(mask):
    """Lv2/Lv7: number of threshold-exceeding connected components (EMERGENCE_LEVELS.md Level 2
    "局在構造の数（閾値超の連結成分数）" / Level 7 droplet/offspring counting). 2D/3D via
    scipy.ndimage.label with full (edge+corner) connectivity.

    mask: boolean array. Returns (count, labeled_array, sizes) where sizes[i] is the pixel/voxel
    count of component i+1 (1-indexed labels), sizes == [] if count == 0.
    """
    m = np.asarray(mask, dtype=bool)
    structure = np.ones((3,) * m.ndim, dtype=int)
    labeled, count = ndimage.label(m, structure=structure)
    sizes = ndimage.sum(m, labeled, index=np.arange(1, count + 1)) if count else np.array([])
    return int(count), labeled, sizes


def droplet_count_and_size(field, thr=0.0):
    """Lv7: droplet count and mean droplet size (EMERGENCE_LEVELS.md Level 7 growth/division:
    counting discrete droplets/offspring by thresholding a scalar order-parameter field, e.g.
    Cahn-Hilliard phi > thr for the high-phi phase).

    Returns (count, mean_size).
    """
    mask = np.asarray(field) > thr
    count, _, sizes = connected_components(mask)
    mean_size = float(sizes.mean()) if count else 0.0
    return count, mean_size


def kinetic_energy(u_components):
    """Lv3: mean kinetic energy 0.5*<|u|^2> (EMERGENCE_LEVELS.md Level 3 "自発運動").
    u_components: sequence of velocity-component arrays [u_x, u_y, (u_z)], all the same shape.
    """
    comps = [np.asarray(c) for c in u_components]
    speed_sq = sum(c ** 2 for c in comps)
    return float(0.5 * np.mean(speed_sq))


def circulation(u_components):
    """Lv3: circulation / vorticity proxy (EMERGENCE_LEVELS.md Level 3 "循環積分 ∮v·dl、渦度").
    2D velocity (u_x, u_y): returns mean |vorticity| = |d(u_y)/dx - d(u_x)/dy| (periodic finite
    difference, unit spacing). 3D velocity (u_x, u_y, u_z): returns mean |curl u|.
    """
    comps = [np.asarray(c) for c in u_components]

    def d(f, axis):
        return np.roll(f, -1, axis) - f

    if len(comps) == 2:
        ux, uy = comps
        vorticity = d(uy, 0) - d(ux, 1)
        return float(np.mean(np.abs(vorticity)))
    if len(comps) == 3:
        ux, uy, uz = comps
        curl_x = d(uz, 1) - d(uy, 2)
        curl_y = d(ux, 2) - d(uz, 0)
        curl_z = d(uy, 0) - d(ux, 1)
        mag = np.sqrt(curl_x ** 2 + curl_y ** 2 + curl_z ** 2)
        return float(np.mean(mag))
    raise ValueError("circulation expects 2 (2D) or 3 (3D) velocity components, got %d" % len(comps))


def occupied_fraction(mask):
    """Lv6: occupied fraction of a boolean structure mask (EMERGENCE_LEVELS.md Level 6 load-bearing
    closure test: compare occupied_fraction between an intact run and an ablation control -- "one
    relation broken collapses the whole" shows up as occupied_fraction collapsing toward 0).
    """
    m = np.asarray(mask, dtype=bool)
    return float(np.mean(m)) if m.size else 0.0
