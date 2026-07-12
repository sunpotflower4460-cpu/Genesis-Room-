#!/usr/bin/env python3
"""渦の共生・橋・三角形 の共通基盤: 3D Gross-Pitaevskii (GPE) 渦線ソルバ。

    i ∂ψ/∂t = [ -½∇² + g|ψ|² - μ ] ψ            (保存的・散逸なし)
    (i - γ) ∂ψ/∂t = [ -½∇² + g|ψ|² - μ ] ψ       (散逸あり、Pitaevskii現象論的減衰)

減衰版を整理すると ∂ψ/∂t = -β (H-μ)ψ,  β = (γ + i)/(1+γ²),  H = -½∇² + g|ψ|²。
  γ=0 → β=i → ユニタリ(ノルム保存)な保存的GPE。
  γ>0 → 実部 -γ/(1+γ²)(H-μ) が加わり、基底状態へ緩和(渦芯を埋める=対消滅)。

数値法: Strang分割のsplit-step Fourier(運動項はFourier空間で対角、相互作用項は実空間で対角)。
周期境界、単位でない格子間隔dxはk_gridのspacingで扱う。

渦のimprint(規律): 完成形(三角形/橋/Hopfion)は初期条件に置かない。渦の「位置と符号」だけを
位相 exp(i·s·atan2) として入れ、共生・橋・面・次元は出力として測る。背景は一様基底状態
ψ=√(μ/g)(周期一様は厳密に定常)。芯は tanh プロファイルで軽く掘り、音波の発生を抑える。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

MODEL_ID = "gpe_vortex_3d"
DEFAULTS = {"g": 1.0, "mu": 1.0, "gamma": 0.0, "dx": 0.5}


def background_amplitude(p):
    """一様基底状態の振幅 √(μ/g)。周期一様は H-μ=0 で厳密に定常。"""
    return np.sqrt(p["mu"] / p["g"])


def make_ground_state(shape, p):
    psi0 = background_amplitude(p)
    return np.full(shape, psi0, dtype=np.complex128)


def imprint_straight_line(psi, x0, y0, charge, p, axis_pair=(2, 1), core_dig=True):
    """z軸に平行な直線渦を(x0,y0)に符号chargeでimprint(位相のみ+芯を軽く掘る)。
    axis_pair=(ax_x, ax_y): 位相を張る2軸(既定 x=axis2, y=axis1、z=axis0が渦線方向)。"""
    ax_x, ax_y = axis_pair
    dx = p["dx"]
    coords = np.indices(psi.shape).astype(float) * dx
    xx = coords[ax_x]
    yy = coords[ax_y]
    theta = np.arctan2(yy - y0, xx - x0)
    psi = psi * np.exp(1j * charge * theta)
    if core_dig:
        r = np.sqrt((xx - x0) ** 2 + (yy - y0) ** 2)
        xi = 1.0 / np.sqrt(2.0 * p["mu"])  # 概略の healing length
        psi = psi * np.tanh(r / xi)
    return psi


def imprint_lines(shape, specs, p, core_dig=True):
    """複数の直線渦をimprint。specs: [(x0,y0,charge), ...]。位置と符号だけを入れる。"""
    psi = make_ground_state(shape, p)
    for (x0, y0, charge) in specs:
        psi = imprint_straight_line(psi, x0, y0, charge, p, core_dig=core_dig)
    return psi


def _step_factors(shape, dt, p):
    _, k2 = k_grid(shape, spacing=p["dx"])
    beta = (p["gamma"] + 1j) / (1.0 + p["gamma"] ** 2)
    kin_half = np.exp(-beta * (k2 / 2.0) * (dt / 2.0))  # 半ステップ運動項(Strang)
    return k2, beta, kin_half


def step(psi, dt, p, kin_half, beta):
    # Strang: 運動半 -> 相互作用全 -> 運動半
    psi_hat = np.fft.fftn(psi)
    psi = np.fft.ifftn(kin_half * psi_hat)
    dens = np.abs(psi) ** 2
    psi = psi * np.exp(-beta * (p["g"] * dens - p["mu"]) * dt)
    psi_hat = np.fft.fftn(psi)
    psi = np.fft.ifftn(kin_half * psi_hat)
    return psi


def total_norm(psi, p):
    return float(np.sum(np.abs(psi) ** 2) * p["dx"] ** 3)


def energy(psi, p):
    """GPEエネルギー E = ∫[ ½|∇ψ|² + ½g|ψ|⁴ - μ|ψ|² ] dV (保存量、γ=0で保存)。"""
    _, k2 = k_grid(psi.shape, spacing=p["dx"])
    psi_hat = np.fft.fftn(psi)
    grad_term = 0.5 * np.sum(k2 * np.abs(psi_hat) ** 2) / psi.size  # Parseval(正規化)
    dens = np.abs(psi) ** 2
    inter = 0.5 * p["g"] * np.sum(dens ** 2) - p["mu"] * np.sum(dens)
    return float((grad_term + inter) * p["dx"] ** 3)


# --------------------------------------------------------------------------
# 渦検出: z断面での平面プラケット巻き数から芯の位置と符号を得る
# --------------------------------------------------------------------------

def detect_vortices_slice(psi_slice, p, amp_frac=0.0):
    """2Dスライス psi_slice(ny,nx)の位相特異点(点渦)を検出。
    返り値: [(x, y, charge), ...] 物理座標(dx込み)、chargeは整数巻き。

    GPEの背景は(乱れた領域でなく)コヒーレントに秩序化しており、位相は芯以外で well-defined。
    振幅マスクは disorderノイズ抑制用で、ここでは逆に「芯そのもの(密度ディップ)」を消してしまう
    ため amp_frac=0.0(実質マスクなし。位相巻き±2πだけで検出)を既定にする。"""
    theta = np.angle(psi_slice)
    amp = np.abs(psi_slice)
    thr = amp_frac * float(amp.max()) if amp.size else 0.0
    w = diag._plaquette_winding(theta, 0, 1, amp, thr)  # (ny,nx)平面
    ys, xs = np.nonzero(w)
    out = []
    for yi, xi in zip(ys, xs):
        # プラケット(yi,xi)は角 (yi,xi)-(yi,xi+1)-(...)-() の中心 ~ (yi+0.5, xi+0.5)
        out.append(((xi + 0.5) * p["dx"], (yi + 0.5) * p["dx"], int(w[yi, xi])))
    return out


def cluster_cores(cores, merge_dist):
    """近接する同符号プラケットを1つの芯にまとめる。cores:[(x,y,q)]。
    返り値: [(x_mean, y_mean, q_sum, n_plaq)]。"""
    remaining = list(cores)
    clusters = []
    while remaining:
        x0, y0, q0 = remaining.pop(0)
        group = [(x0, y0, q0)]
        changed = True
        while changed:
            changed = False
            keep = []
            for (x, y, q) in remaining:
                if any(np.hypot(x - gx, y - gy) <= merge_dist and np.sign(q) == np.sign(gq)
                       for (gx, gy, gq) in group):
                    group.append((x, y, q))
                    changed = True
                else:
                    keep.append((x, y, q))
            remaining = keep
        xs = np.mean([g[0] for g in group])
        ys = np.mean([g[1] for g in group])
        qsum = sum(g[2] for g in group)
        clusters.append((float(xs), float(ys), int(np.sign(qsum)) if qsum != 0 else 0, len(group)))
    return clusters


def slice_core_summary(psi, p, z_index=None, merge_dist=None):
    """中央(既定)z断面での芯: 個数, 正芯数, 負芯数, 位置リスト, 総巻き。"""
    if z_index is None:
        z_index = psi.shape[0] // 2
    merge_dist = merge_dist or 3.0 * p["dx"]
    raw = detect_vortices_slice(psi[z_index], p)
    clusters = cluster_cores(raw, merge_dist)
    n_pos = sum(1 for c in clusters if c[2] > 0)
    n_neg = sum(1 for c in clusters if c[2] < 0)
    net_winding = sum(c[2] for c in clusters)
    return {"n_cores": len(clusters), "n_pos": n_pos, "n_neg": n_neg,
            "net_winding": net_winding, "cores": clusters}


def line_length_proxy(psi, amp_frac=0.0):
    """3D全体の渦線長プロキシ(全3面配向のプラケット巻きの非零数)。
    GPEでは芯を消さないよう amp_frac=0.0(位相巻きのみ)を既定にする。"""
    return diag.winding_defect_count(psi, amplitude_threshold_frac=amp_frac)


def run(shape, specs, steps, dt, p=None, snapshot_every=None, core_dig=True,
        z_index=None, merge_dist=None):
    """直線渦configをimprintしてGPEで発展。芯の時系列を返す。"""
    params = dict(DEFAULTS)
    if p:
        params.update(p)
    psi = imprint_lines(shape, specs, params, core_dig=core_dig)
    _, beta, kin_half = _step_factors(shape, dt, params)
    snapshot_every = snapshot_every or max(1, steps // 40)
    norm0 = total_norm(psi, params)
    e0 = energy(psi, params)
    snapshots = []
    diverged = False
    for t in range(steps):
        psi = step(psi, dt, params, kin_half, beta)
        if not np.all(np.isfinite(psi)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            summ = slice_core_summary(psi, params, z_index=z_index, merge_dist=merge_dist)
            summ["step"] = t
            summ["t"] = round(t * dt, 3)
            summ["line_len"] = line_length_proxy(psi)
            summ["norm"] = total_norm(psi, params)
            summ["energy"] = energy(psi, params)
            summ["min_density"] = float((np.abs(psi) ** 2).min())
            snapshots.append(summ)
    phys = {"norm0": norm0, "energy0": e0, "diverged": diverged,
            "norm_final": total_norm(psi, params), "energy_final": energy(psi, params)}
    return snapshots, phys, psi


if __name__ == "__main__":
    # 簡易セルフテスト: 単一渦は定常のはず、逆符号ペアは並進、減衰で対消滅。
    p = dict(DEFAULTS, dx=0.5)
    shape = (48, 48, 48)
    box = shape[1] * p["dx"]
    print("box size =", box, "healing xi ~", 1 / np.sqrt(2 * p["mu"]))
    snaps, phys, psi = run(shape, [(box * 0.4, box * 0.5, +1), (box * 0.6, box * 0.5, -1)],
                           steps=200, dt=0.02, p=p)
    for s in snaps[::5]:
        print("t=%.2f n_cores=%d (+%d/-%d) net=%d line=%d norm=%.3f E=%.3f"
              % (s["t"], s["n_cores"], s["n_pos"], s["n_neg"], s["net_winding"],
                 s["line_len"], s["norm"], s["energy"]))
