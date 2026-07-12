#!/usr/bin/env python3
"""genesis/ic_families_real.py -- 実数場(Gray-Scott の v・Model H の phi)向けの7 family。

genesis/ic_families_3d.py（複素場CGL向け）と同じ7つの切り口を、実数場に翻訳する。
「位相」は実数場には存在しないため、複素場での「位相多様性」は実数場では「符号多様性」に
翻訳する（同じ役割: 各構造が独立にバラバラの"向き"を持つか、全て同じ"向き"かを対比する）。
第8監査は複素版と同じ: 完成形を置かない、種・ノイズ・スペクトル・幾何構造のみ。

各生成器は shape と numpy Generator を受け取り、平均0付近の実数ND配列（摂動場、まだ
"背景+この摂動"には合成していない）を返す——背景への合成は呼び出し側(law adapter)が行う。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from genesis.solvers import k_grid  # noqa: E402

FAMILIES = ["spectral_powerlaw", "bandpass", "white_amp", "seeds_signed", "line_or_ring_patches",
            "patches", "seeds"]

PARAM_RANGES = {
    "spectral_powerlaw": {"beta": (0.5, 4.0), "amp": (0.1, 0.6)},
    "bandpass": {"klo": (0.3, 4.0), "kw": (0.5, 3.0), "amp": (0.1, 0.6)},
    "white_amp": {"amp": (0.1, 0.6)},
    "seeds_signed": {"n_seeds": (2, 8), "width": (2.0, 6.0), "amp": (0.2, 0.6)},
    "line_or_ring_patches": {"n_defects": (1, 3), "amp": (0.2, 0.6)},
    "patches": {"n_patches": (1, 6), "radius_lo": (2.5, 6.0), "radius_hi": (6.0, 12.0),
                "amp": (0.2, 0.6)},
    "seeds": {"n_seeds": (2, 8), "width": (2.0, 6.0), "amp": (0.2, 0.6)},
}


def _wrapped_dist2(shape, center):
    grids = np.indices(shape).astype(float)
    d2 = np.zeros(shape)
    for a in range(len(shape)):
        d = np.abs(grids[a] - center[a])
        d = np.minimum(d, shape[a] - d)
        d2 += d ** 2
    return d2


def sample_params(family, rng):
    ranges = PARAM_RANGES[family]
    p = {}
    for key, (lo, hi) in ranges.items():
        if key in ("n_seeds", "n_defects", "n_patches"):
            p[key] = int(rng.integers(lo, hi + 1))
        else:
            p[key] = float(rng.uniform(lo, hi))
    return p


def spectral_powerlaw(shape, rng, beta=2.0, amp=0.3):
    """べき則スペクトルの実数ノイズ P(k)~|k|^-beta。相関長は法則側の空間スケール(健全な
    パターン間隔)に近いほど点火しやすい(GSで予備確認済み)。"""
    noise = rng.standard_normal(shape)
    noise_hat = np.fft.fftn(noise)
    _, k2 = k_grid(shape)
    k = np.sqrt(k2)
    filt = np.zeros_like(k)
    mask = k > 0
    filt[mask] = k[mask] ** (-beta / 2.0)
    field = np.real(np.fft.ifftn(noise_hat * filt))
    field = field - field.mean()
    scale = np.std(field) + 1e-12
    return field / scale * amp


def bandpass(shape, rng, klo=1.0, kw=2.0, amp=0.3):
    noise = rng.standard_normal(shape)
    noise_hat = np.fft.fftn(noise)
    _, k2 = k_grid(shape)
    k = np.sqrt(k2)
    mask = (k >= klo) & (k < klo + kw)
    field = np.real(np.fft.ifftn(noise_hat * mask))
    field = field - field.mean()
    scale = np.std(field) + 1e-12
    return field / scale * amp


def white_amp(shape, rng, amp=0.3):
    """空間相関のない白色ノイズ(GSでは点火しにくいと予想、決定的な家系間対比)。"""
    return amp * rng.standard_normal(shape)


def seeds_signed(shape, rng, n_seeds=5, width=3.0, amp=0.4):
    """符号付き種: 各バンプが独立なランダム符号(+/-)を持つ(複素版のseeds_phaseの実数翻訳、
    "位相多様性"->"符号多様性")。"""
    field = np.zeros(shape)
    for _ in range(n_seeds):
        center = [rng.uniform(0, s) for s in shape]
        d2 = _wrapped_dist2(shape, center)
        bump = amp * np.exp(-d2 / (2 * width ** 2))
        sign = rng.choice([-1.0, 1.0])
        field = field + sign * bump
    return field


def line_or_ring_patches(shape, rng, n_defects=1, amp=0.4):
    """3D特有の幾何構造: 渦線/渦リングと同じ形状(周期箱を貫く管 or 局所トーラス)だが、
    実数場には位相がないため巻き数は持たない(振幅の帯のみ)。genesis/ic_families_3d.py の
    vortex_charges から位相項を落とした幾何学的アナログ。"""
    field = np.zeros(shape)
    grids = np.indices(shape).astype(float)
    for _ in range(n_defects):
        mode = rng.choice(["line", "loop"])
        if mode == "line":
            axis = int(rng.integers(0, len(shape)))
            other = [a for a in range(len(shape)) if a != axis]
            x0 = rng.uniform(0, shape[other[0]])
            y0 = rng.uniform(0, shape[other[1]])
            d0 = np.abs(grids[other[0]] - x0)
            d0 = np.minimum(d0, shape[other[0]] - d0)
            d1 = np.abs(grids[other[1]] - y0)
            d1 = np.minimum(d1, shape[other[1]] - d1)
            r = np.sqrt(d0 ** 2 + d1 ** 2)
            field = field + amp * np.exp(-(r / 2.5) ** 2)
        else:
            axis = int(rng.integers(0, len(shape)))
            other = [a for a in range(len(shape)) if a != axis]
            center = [rng.uniform(0, s) for s in shape]
            R = rng.uniform(3.0, max(min(shape) / 4.0, 4.0))
            r_tube = rng.uniform(1.5, 3.0)
            d_axis = np.abs(grids[axis] - center[axis])
            d_axis = np.minimum(d_axis, shape[axis] - d_axis)
            d0 = np.abs(grids[other[0]] - center[other[0]])
            d0 = np.minimum(d0, shape[other[0]] - d0)
            d1 = np.abs(grids[other[1]] - center[other[1]])
            d1 = np.minimum(d1, shape[other[1]] - d1)
            rho = np.sqrt(d0 ** 2 + d1 ** 2)
            dist_to_tube = np.sqrt((rho - R) ** 2 + d_axis ** 2)
            field = field + amp * np.exp(-(dist_to_tube / r_tube) ** 2)
    return field


def patches(shape, rng, n_patches=3, radius_lo=3.0, radius_hi=8.0, amp=0.4):
    """滑らかなドメインのパッチ(全て同じ符号、パッチ間で独立でない)。"""
    field = np.zeros(shape)
    for _ in range(n_patches):
        center = [rng.uniform(0, s) for s in shape]
        radius = rng.uniform(radius_lo, radius_hi)
        d2 = _wrapped_dist2(shape, center)
        field = field + amp * np.exp(-d2 / (2 * radius ** 2))
    return field


def seeds(shape, rng, n_seeds=5, width=3.0, amp=0.4):
    """符号なし種(対照・全て正、多様性なし)。"""
    field = np.zeros(shape)
    for _ in range(n_seeds):
        center = [rng.uniform(0, s) for s in shape]
        d2 = _wrapped_dist2(shape, center)
        field = field + amp * np.exp(-d2 / (2 * width ** 2))
    return field


GENERATORS = {
    "spectral_powerlaw": spectral_powerlaw, "bandpass": bandpass, "white_amp": white_amp,
    "seeds_signed": seeds_signed, "line_or_ring_patches": line_or_ring_patches,
    "patches": patches, "seeds": seeds,
}


def make_ic(family, shape, rng, params=None):
    p = params if params is not None else sample_params(family, rng)
    return GENERATORS[family](shape, rng, **p), p


if __name__ == "__main__":
    shape = (24, 24, 24)
    rng = np.random.default_rng(0)
    for fam in FAMILIES:
        field, p = make_ic(fam, shape, rng)
        print("%-20s params=%s  mean=%.4f std=%.4f max=%.3f"
              % (fam, {k: round(v, 3) if isinstance(v, float) else v for k, v in p.items()},
                 field.mean(), field.std(), field.max()))
