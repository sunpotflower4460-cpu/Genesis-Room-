#!/usr/bin/env python3
"""依頼A/B/C共通基盤: 区画化(compartmentalized)genesis — 生産を液滴"内部"でのみ行う二場モデル(3D)。

    du/dt = M*lap(mu) + R,   mu = f'(u) - kappa*lap(u),  f'(u)=2*W*u*(1-u)*(1-2u)
    dv/dt = D_v*lap(v) - R
    R = kp*v*chi(u) - kd*chi(u)^2                       (区画化: 生産は内部のみ)
    R_seeded = kp_seed*v*(1-chi(u)) + kp_in*v*chi(u) - kd*chi(u)^2   (依頼B: 弱い外部seed+強い内部区画化)
    chi(u) = 0.5*(1+tanh((u-0.5)/w))

区画化の物理的機構: 希薄背景(chi≈0)では R≈kp*v*0 - kd*0 ≈ 0（生産も分解もほぼ起きない、
拡散でしか変化しない）。液滴内部(chi≈1)でのみ v→u 変換が起きる。これが「なぜ膜が要るか」
——生産を外部リザーバから隔離することで、背景を汚さずに個体を維持できる。

Round2の知見（chi^2ゲートは"安全"kdで液滴を救ったが独自の安全window再較正が必要、との
発見）を受けて、この区画化モデルでは生産項もchiでゲートし直した——理論上、背景では生産も
分解もほぼ0になるはずで、H-B1的な「安全比」を探す必要自体がなくなる、というのが依頼Aの核心
仮説。
"""

import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import diagnostics as diag  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

MODEL_ID = "compartmentalized_genesis_3d"
DEFAULTS = {"M": 1.0, "kappa": 1.0, "W": 1.0, "w": 0.15, "D_v": 2.0, "k_p": 0.05, "k_d": 0.03}
U_SP_LOW = (3 - np.sqrt(3)) / 6  # ~0.2113, f''(u)=0 (W=1)
U_SP_HIGH = (3 + np.sqrt(3)) / 6  # ~0.7887


def chi(u, w):
    return 0.5 * (1.0 + np.tanh((u - 0.5) / w))


def reaction_compartmentalized(u, v, p):
    """依頼A: R = kp*v*chi(u) - kd*chi(u)^2 -- 生産は内部のみ。"""
    c = chi(u, p["w"])
    return p["k_p"] * v * c - p["k_d"] * c ** 2


def reaction_seeded_compartmentalized(u, v, p):
    """依頼B: 弱い外部seed(1-chi) + 強い内部区画化(chi)。p に k_p_seed, k_p_in が必要。"""
    c = chi(u, p["w"])
    return p["k_p_seed"] * v * (1.0 - c) + p["k_p_in"] * v * c - p["k_d"] * c ** 2


def make_uniform_initial(shape, u_mean, rng, u_noise=0.02, v_mean=0.5):
    """一様(u_mean)+微小ノイズ。液滴は入れない（NO seeded droplets, 依頼A/Bの基本形）。"""
    u = u_mean + u_noise * rng.standard_normal(shape)
    v = np.full(shape, v_mean, dtype=float)
    return u, v


def step(u, v, dt, p, k2, reaction_fn=reaction_compartmentalized):
    R = reaction_fn(u, v, p)
    fprime = 2.0 * p["W"] * u * (1.0 - u) * (1.0 - 2.0 * u)
    fprime_hat = np.fft.fftn(fprime)
    R_hat = np.fft.fftn(R)
    uhat = np.fft.fftn(u)
    vhat = np.fft.fftn(v)
    denom_u = 1.0 + dt * p["M"] * p["kappa"] * k2 ** 2
    uhat_new = (uhat + dt * (-p["M"] * k2 * fprime_hat + R_hat)) / denom_u
    denom_v = 1.0 + dt * p["D_v"] * k2
    vhat_new = (vhat - dt * R_hat) / denom_v
    return np.real(np.fft.ifftn(uhat_new)), np.real(np.fft.ifftn(vhat_new))


def _run_core(u0, v0, shape, steps, dt, params, snapshot_every, reaction_fn):
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    _, k2 = k_grid(shape)
    u, v = u0, v0
    mass0 = float(np.mean(u + v))
    snapshots = []
    diverged = False
    u_bg_series = []
    for t in range(steps):
        u, v = step(u, v, dt, p, k2, reaction_fn=reaction_fn)
        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
            diverged = True
            break
        dilute = u < 0.5
        u_bg_series.append(float(np.mean(u[dilute])) if dilute.any() else float(np.mean(u)))
        if t % snapshot_every == 0 or t == steps - 1:
            n, _, sizes = diag.connected_components(u > 0.5)
            snapshots.append({"step": t, "u": u.copy(), "v": v.copy(), "n_droplets": n,
                              "sizes": sorted(sizes.tolist(), reverse=True)[:10]})
    mass1 = float(np.mean(u + v))
    if not snapshots:
        n, _, sizes = diag.connected_components(u > 0.5)
        snapshots = [{"step": 0, "u": u.copy(), "v": v.copy(), "n_droplets": n,
                     "sizes": sorted(sizes.tolist(), reverse=True)[:10]}]
    phys = {"mass_drift": mass1 - mass0, "diverged": diverged, "u_bg_series": u_bg_series,
            "u_bg_final": u_bg_series[-1] if u_bg_series else None,
            "u_bg_max": max(u_bg_series) if u_bg_series else None}
    return snapshots, phys


def run(shape, steps, dt, seed, params=None, snapshot_every=None, u_mean=0.35, v_mean=0.5,
        u_noise=0.02, reaction_fn=reaction_compartmentalized):
    """natural genesis: t=0（一様u_mean+微小ノイズ、完成した液滴なし）から中断なく発展させる。"""
    rng = np.random.default_rng(seed)
    u0, v0 = make_uniform_initial(shape, u_mean, rng, u_noise=u_noise, v_mean=v_mean)
    snapshot_every = snapshot_every or max(1, steps // 30)
    return _run_core(u0, v0, shape, steps, dt, params, snapshot_every, reaction_fn)


def run_droplet_probe(shape, R0, steps, dt, seed, params=None, snapshot_every=None,
                       u_background=0.15, v_uniform=0.5, reaction_fn=reaction_compartmentalized,
                       dx=1.0):
    """決定的対照/診断: 半径R0の単一液滴（依頼Cの分裂テスト用）。"""
    from genesis.dividing_protocell_3d import make_droplet_initial
    rng = np.random.default_rng(seed)
    u0, v0 = make_droplet_initial(shape, R0, dx=dx, u_background=u_background, v_uniform=v_uniform,
                                   rng=rng)
    snapshot_every = snapshot_every or max(1, steps // 30)
    return _run_core(u0, v0, shape, steps, dt, params, snapshot_every, reaction_fn)
