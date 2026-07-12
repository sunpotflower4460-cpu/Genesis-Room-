#!/usr/bin/env python3
"""genesis/ic_families_3d.py -- IC family を大量ランダム生成する（法則は変えない、変えるのは
初期条件だけ）。指示書「大量ランダムIC探索」の7 family を3D化。

第8監査（各 family が守ること）: 完成形（安定した渦線構造・分裂後の2つ・木構造など）を
初期条件に置かない。置くのは「種・ノイズ・対称構造・位相」だけ——どこまで自然に育つかは
run が出力として測る。振幅は一様近ゼロ+ノイズのオーダー（noise_floor 程度〜amp 数割）に
抑え、法則(G001 CGL)の定常振幅|A|~1に対して常に「小さい種」であることを保証する。

全 family は shape(3D) と numpy Generator を受け取り、同じ shape の complex128 配列を返す。
乱数はすべて渡された rng からのみ引くため、同じ rng 状態列（=同じ seed）なら常に同じ IC が
決定的に再現される。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from genesis.solvers import k_grid  # noqa: E402

FAMILIES = ["spectral_powerlaw", "bandpass", "white_amp", "seeds_phase", "vortex_charges",
            "patches", "seeds"]

# 第8監査: ここに書くのは「IC生成パラメータのランダム化レンジ」であって完成形ではない。
# search_space.yaml にも同じレンジを明文化してある。
PARAM_RANGES = {
    "spectral_powerlaw": {"beta": (0.5, 4.0), "amp": (0.01, 0.3)},
    "bandpass": {"klo": (0.5, 8.0), "kw": (0.5, 4.0), "amp": (0.01, 0.3)},
    "white_amp": {"amp": (0.005, 0.5)},
    "seeds_phase": {"n_seeds": (2, 10), "width": (1.5, 6.0), "amp": (0.1, 0.8)},
    "vortex_charges": {"n_defects": (1, 4), "amp": (0.1, 0.6)},
    "patches": {"n_patches": (2, 8), "radius_lo": (2.0, 6.0), "radius_hi": (6.0, 14.0),
                "amp": (0.1, 0.6)},
    "seeds": {"n_seeds": (2, 10), "width": (1.5, 6.0), "amp": (0.1, 0.8)},
}

NOISE_FLOOR = 0.01  # 全familyに共通の背景ノイズ床(一様近ゼロ+ノイズという法則側の前提を保つ)


def _wrapped_dist2(shape, center):
    """周期境界での最小像距離の二乗（3Dボックスの端をまたぐ種も正しく扱う）。"""
    grids = np.indices(shape).astype(float)
    d2 = np.zeros(shape)
    for a in range(len(shape)):
        d = np.abs(grids[a] - center[a])
        d = np.minimum(d, shape[a] - d)
        d2 += d ** 2
    return d2


def sample_params(family, rng):
    """family の PARAM_RANGES からランダムにパラメータを1組引く。"""
    ranges = PARAM_RANGES[family]
    p = {}
    for key, (lo, hi) in ranges.items():
        if key in ("n_seeds", "n_defects", "n_patches"):
            p[key] = int(rng.integers(lo, hi + 1))
        else:
            p[key] = float(rng.uniform(lo, hi))
    return p


def spectral_powerlaw(shape, rng, beta=2.0, amp=0.05):
    """べき則スペクトルのノイズ: P(k) ~ |k|^-beta。scale-free的な始まり(予備調査で上位)。"""
    noise = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    noise_hat = np.fft.fftn(noise)
    _, k2 = k_grid(shape)
    k = np.sqrt(k2)
    filt = np.zeros_like(k)
    mask = k > 0
    filt[mask] = k[mask] ** (-beta / 2.0)
    field = np.fft.ifftn(noise_hat * filt)
    field = field - field.mean()
    scale = np.std(np.abs(field)) + 1e-12
    return (field / scale * amp).astype(np.complex128)


def bandpass(shape, rng, klo=2.0, kw=2.0, amp=0.05):
    """白色ノイズを波数帯[klo, klo+kw]だけ通す(初期の長さスケールを制御、予備調査で上位)。"""
    noise = rng.standard_normal(shape) + 1j * rng.standard_normal(shape)
    noise_hat = np.fft.fftn(noise)
    _, k2 = k_grid(shape)
    k = np.sqrt(k2)
    mask = (k >= klo) & (k < klo + kw)
    field = np.fft.ifftn(noise_hat * mask)
    field = field - field.mean()
    scale = np.std(np.abs(field)) + 1e-12
    return (field / scale * amp).astype(np.complex128)


def white_amp(shape, rng, amp=0.05):
    """振幅ランダムの白色ノイズ(既存 make_initial の amp をランダム化した版)。"""
    return (amp * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))).astype(np.complex128)


def seeds_phase(shape, rng, n_seeds=5, width=3.0, amp=0.4):
    """位相付き種: 各バンプが独立なランダム位相を持つ(予備調査: 位相なしより良い)。"""
    field = NOISE_FLOOR * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
    for _ in range(n_seeds):
        center = [rng.uniform(0, s) for s in shape]
        d2 = _wrapped_dist2(shape, center)
        bump = amp * np.exp(-d2 / (2 * width ** 2))
        phase = rng.uniform(0, 2 * np.pi)
        field = field + bump * np.exp(1j * phase)
    return field.astype(np.complex128)


def vortex_charges(shape, rng, n_defects=1, amp=0.3, charge_choices=(-2, -1, 1, 2)):
    """トポロジー的IC: 3D渦線(周期箱を貫く直線=torus上のループと等価)または渦リング(局所ループ)
    を位相にimprintする(渦"を撒く"だけ。「橋/三角形にせよ」等の完成形は置かない)。"""
    field = NOISE_FLOOR * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
    grids = np.indices(shape).astype(float)
    for _ in range(n_defects):
        mode = rng.choice(["line", "loop"])
        charge = int(rng.choice(charge_choices))
        if mode == "line":
            axis = int(rng.integers(0, len(shape)))
            other = [a for a in range(len(shape)) if a != axis]
            x0 = rng.uniform(0, shape[other[0]])
            y0 = rng.uniform(0, shape[other[1]])
            theta = np.arctan2(grids[other[1]] - y0, grids[other[0]] - x0)
            phase = charge * theta
            field = field + amp * np.exp(1j * phase)
        else:  # loop: 局所的な渦リング(トーラス状の管の中でだけ位相が巻く)
            axis = int(rng.integers(0, len(shape)))
            other = [a for a in range(len(shape)) if a != axis]
            center = [rng.uniform(0, s) for s in shape]
            R = rng.uniform(3.0, max(min(shape) / 4.0, 4.0))
            r_tube = rng.uniform(1.5, 3.5)
            d_axis = np.abs(grids[axis] - center[axis])
            d_axis = np.minimum(d_axis, shape[axis] - d_axis)
            d0 = np.abs(grids[other[0]] - center[other[0]])
            d0 = np.minimum(d0, shape[other[0]] - d0)
            d1 = np.abs(grids[other[1]] - center[other[1]])
            d1 = np.minimum(d1, shape[other[1]] - d1)
            rho = np.sqrt(d0 ** 2 + d1 ** 2)
            dist_to_tube = np.sqrt((rho - R) ** 2 + d_axis ** 2)
            phase = charge * np.arctan2(d_axis, rho - R)
            envelope = amp * np.exp(-np.maximum(dist_to_tube - r_tube, 0) ** 2 / (2 * 2.0 ** 2))
            field = field + envelope * np.exp(1j * phase)
    return field.astype(np.complex128)


def patches(shape, rng, n_patches=4, radius_lo=4.0, radius_hi=10.0, amp=0.3):
    """滑らかなドメインのパッチ: 振幅の塊、位相は全体で1つ(パッチ間で独立でない)。"""
    global_phase = rng.uniform(0, 2 * np.pi)
    field = NOISE_FLOOR * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
    for _ in range(n_patches):
        center = [rng.uniform(0, s) for s in shape]
        radius = rng.uniform(radius_lo, radius_hi)
        d2 = _wrapped_dist2(shape, center)
        bump = amp * np.exp(-d2 / (2 * radius ** 2))
        field = field + bump * np.exp(1j * global_phase)
    return field.astype(np.complex128)


def seeds(shape, rng, n_seeds=5, width=3.0, amp=0.4):
    """位相なし種(対照・予備調査で最悪=0%): 実正の振幅バンプのみ、位相多様性なし。"""
    field = NOISE_FLOOR * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))
    for _ in range(n_seeds):
        center = [rng.uniform(0, s) for s in shape]
        d2 = _wrapped_dist2(shape, center)
        bump = amp * np.exp(-d2 / (2 * width ** 2))
        field = field + bump  # phase=0, no diversity by construction
    return field.astype(np.complex128)


GENERATORS = {
    "spectral_powerlaw": spectral_powerlaw, "bandpass": bandpass, "white_amp": white_amp,
    "seeds_phase": seeds_phase, "vortex_charges": vortex_charges, "patches": patches,
    "seeds": seeds,
}


def make_ic(family, shape, rng, params=None):
    """family名から IC を生成。params が None ならその場でランダムサンプル(sample_paramsも
    同じ rng から引くので、rng の状態列だけで完全に再現可能)。"""
    p = params if params is not None else sample_params(family, rng)
    return GENERATORS[family](shape, rng, **p), p


if __name__ == "__main__":
    shape = (24, 24, 24)
    rng = np.random.default_rng(0)
    for fam in FAMILIES:
        field, p = make_ic(fam, shape, rng)
        print("%-18s params=%s  |A|mean=%.4f |A|max=%.4f mean(A)=%.2e"
              % (fam, {k: round(v, 3) if isinstance(v, float) else v for k, v in p.items()},
                 np.abs(field).mean(), np.abs(field).max(), abs(field.mean())))
