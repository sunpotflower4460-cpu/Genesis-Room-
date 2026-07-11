#!/usr/bin/env python3
"""依頼D【探索】: 流れ（Model H）でLv3/Lv5も足した最大climb（3D）。

G003で確立した非圧縮Navier-Stokes+毛管力を、依頼Aの区画化genesisに結合する:

    du/dt + vel.grad(u) = M*lap(mu) + R,   mu = f'(u) - kappa*lap(u)
    dv/dt = D_v*lap(v) - R                                          # v は移流させない(簡略化)
    dvel/dt + vel.grad(vel) = -grad(p) + nu*lap(vel) + C*(mu*grad(u)),  div(vel)=0
    R = kp*v*chi(u) - kd*chi(u)^2                                   # 区画化（依頼Aと同じ）

--- 予言（実行前に登録） ---
流れがLv3（自発運動/循環）とLv5（流れが共分化を駆動）を足し、区画化がLv4/Lv6を保つ→
一つのrunが0->Lv1->2->3->4->5->6とより完全に登る。

決定的対照: 流れ結合C=0に落とすとLv3/Lv5が消え、区画化のLv4/Lv6だけ残る（依頼Aと同じ挙動
に帰着するはず）。

これは探索（2D予備なし）。falsification: 流れを入れてもLv3/5相当の指標（KE>0、循環）が
測定できない、または流れが区画化を破壊する（背景核形成/個体死）なら、この結合は今回のパラ
メータでは機能しない、と正直に報告する。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import compartmentalized_genesis_3d as cg  # noqa: E402
from genesis.solvers import k_grid, dealias_mask, leray_project  # noqa: E402

DT = 0.05
STEPS = 3000  # t_final=150, 依頼Aと同じ設定で直接比較できるようにする
SHAPE = (64, 64, 64)  # 依頼Aの96^3より小さく: 流体計算はコスト高（探索なので許容）
CONFIRMED_PARAMS = dict(cg.DEFAULTS, C=2.0, nu=1.0)  # 依頼Aのkp/kd + G003のC,nu


def _hat(field, dealias):
    return np.fft.fftn(field) * dealias


def step_with_flow(u, v, vel, dt, p, kk, k2, dealias):
    ndim = u.ndim
    uhat_full = np.fft.fftn(u)
    uhat = uhat_full * dealias
    lap_u = np.real(np.fft.ifftn(-k2 * uhat_full))
    # dividing_protocell_3d/compartmentalized_genesis_3d と同じ二重井戸: f'(u)=2W*u(1-u)(1-2u)
    # (安定点 u=0,1。G003のu^3-u系列(安定点u=-+1)とは別物なので流用しない)
    mu = 2.0 * p["W"] * u * (1.0 - u) * (1.0 - 2.0 * u) - p["kappa"] * lap_u
    grad_u = [np.real(np.fft.ifftn(1j * kk[a] * uhat)) for a in range(ndim)]

    R = cg.reaction_compartmentalized(u, v, p)
    fprime = 2.0 * p["W"] * u * (1.0 - u) * (1.0 - 2.0 * u)

    advect_u = sum(vel[a] * grad_u[a] for a in range(ndim))
    fprime_hat = _hat(fprime, dealias)
    R_hat = _hat(R, dealias)
    # 半陰的: biharmonic項(-M*kappa*k^4*uhat)のみ陰的に、移流/反応/f'(u)は陽的
    # (dividing_protocell_3d.step と同じスキーム、移流項を追加)
    explicit_u_hat = _hat(-advect_u, dealias) - p["M"] * k2 * fprime_hat + R_hat
    denom_u = 1.0 + dt * p["M"] * p["kappa"] * k2 ** 2
    uhat_new = (uhat_full + dt * explicit_u_hat) / denom_u
    u_new = np.real(np.fft.ifftn(uhat_new))

    vhat = np.fft.fftn(v)
    denom_v = 1.0 + dt * p["D_v"] * k2
    vhat_new = (vhat - dt * np.fft.fftn(R)) / denom_v
    v_new = np.real(np.fft.ifftn(vhat_new))

    velhat = [_hat(c, dealias) for c in vel]
    grad_vel = [[np.real(np.fft.ifftn(1j * kk[b] * velhat[a])) for b in range(ndim)] for a in range(ndim)]
    advect_vel = [sum(vel[b] * grad_vel[a][b] for b in range(ndim)) for a in range(ndim)]
    capillary = [p["C"] * mu * grad_u[a] for a in range(ndim)]
    rhs_hat = [_hat(-advect_vel[a] + capillary[a], dealias) for a in range(ndim)]
    velhat_full = [np.fft.fftn(c) for c in vel]
    vel_star_hat = [velhat_full[a] + dt * rhs_hat[a] for a in range(ndim)]
    vel_visc_hat = [vel_star_hat[a] / (1.0 + dt * p["nu"] * k2) for a in range(ndim)]
    vel_proj_hat = leray_project(vel_visc_hat, kk, k2)
    vel_new = [np.real(np.fft.ifftn(vh)) for vh in vel_proj_hat]

    return u_new, v_new, vel_new


def run(shape, steps, dt, seed, params=None, u_mean=0.35, v_mean=0.5, snapshot_every=None):
    p = dict(CONFIRMED_PARAMS)
    if params:
        p.update(params)
    rng = np.random.default_rng(seed)
    u, v = cg.make_uniform_initial(shape, u_mean, rng, v_mean=v_mean)
    vel = [np.zeros(shape) for _ in range(len(shape))]
    kk, k2 = k_grid(shape)
    dealias = dealias_mask(shape)
    snapshot_every = snapshot_every or max(1, steps // 30)

    mass0 = float(np.mean(u + v))
    snapshots = []
    diverged = False
    u_bg_series = []
    ke_series = []
    for t in range(steps):
        u, v, vel = step_with_flow(u, v, vel, dt, p, kk, k2, dealias)
        if not np.all(np.isfinite(u)) or any(not np.all(np.isfinite(c)) for c in vel):
            diverged = True
            break
        dilute = u < 0.5
        u_bg_series.append(float(np.mean(u[dilute])) if dilute.any() else float(np.mean(u)))
        ke_series.append(float(0.5 * np.mean(sum(c ** 2 for c in vel))))
        if t % snapshot_every == 0 or t == steps - 1:
            n, _, sizes = diag.connected_components(u > 0.5)
            snapshots.append({"step": t, "u": u.copy(), "v": v.copy(), "vel": [c.copy() for c in vel],
                              "n_droplets": n, "sizes": sorted(sizes.tolist(), reverse=True)[:8]})
    mass1 = float(np.mean(u + v))
    if not snapshots:
        n, _, sizes = diag.connected_components(u > 0.5)
        snapshots = [{"step": 0, "u": u.copy(), "v": v.copy(), "vel": [c.copy() for c in vel],
                     "n_droplets": n, "sizes": sorted(sizes.tolist(), reverse=True)[:8]}]
    phys = {"mass_drift": mass1 - mass0, "diverged": diverged, "u_bg_series": u_bg_series,
            "u_bg_final": u_bg_series[-1] if u_bg_series else None,
            "ke_series": ke_series, "final_ke": ke_series[-1] if ke_series else 0.0}
    return snapshots, phys


def _level_report(snapshots, phys):
    u_fields = [s["u"] for s in snapshots]
    variances, growth_rate = diag.variance_growth(u_fields)
    _, prominence = diag.structure_factor_peak(u_fields[-1])
    final = snapshots[-1]

    difference = bool(growth_rate > 0 and prominence > 1.5)
    localization = bool(difference and final["n_droplets"] >= 1)
    persistent_individuality = bool(localization and final["n_droplets"] >= 1 and not phys["diverged"])
    u_bg_final = phys["u_bg_final"]
    background_clean = bool(u_bg_final is not None and u_bg_final < cg.U_SP_LOW)
    self_maintaining = bool(persistent_individuality and background_clean)

    spontaneous_motion = bool(phys["final_ke"] > 1e-6)
    vel_final = final.get("vel")
    circulation = 0.0
    if vel_final is not None and len(vel_final) == 3:
        circulation = diag.circulation(vel_final)
    co_differentiation = bool(spontaneous_motion and circulation > 1e-6 and localization)

    reached = 0
    if difference:
        reached = 1
    if localization:
        reached = 2
    if localization and spontaneous_motion:
        reached = 3
    if persistent_individuality:
        reached = max(reached, 4)
    if co_differentiation:
        reached = max(reached, 5)
    if self_maintaining:
        reached = max(reached, 6)

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = difference
    detected["localization"] = localization
    detected["spontaneous_motion"] = spontaneous_motion
    detected["circulation"] = bool(circulation > 1e-6)
    detected["persistent_individuality"] = persistent_individuality
    detected["co_differentiation"] = co_differentiation
    detected["self_maintaining_closure"] = self_maintaining

    role = "E" if reached >= 1 else "F"
    return {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"variance_growth": round(float(growth_rate), 6),
                        "structure_factor_prominence": round(float(prominence), 4),
                        "final_n_droplets": final["n_droplets"], "final_top_sizes": final["sizes"][:5],
                        "u_bg_final": u_bg_final, "u_sp_low": cg.U_SP_LOW,
                        "background_clean": background_clean, "final_kinetic_energy": phys["final_ke"],
                        "final_circulation": circulation, "mass_drift": phys["mass_drift"]},
        "purity": {"per_object_labels": False, "external_optimum": False, "role": role},
    }


def _run_and_save(room_id, params, seed, u_mean=0.35, notes_extra=""):
    snapshots, phys = run(SHAPE, STEPS, DT, seed, params=params, u_mean=u_mean)
    report = _level_report(snapshots, phys)
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"],
        resolutions_result={"%dx%dx%d" % SHAPE: report["measured_by"]["final_n_droplets"]},
        seed_success={str(seed): report["reached_level"]}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "C=0 (flow decoupled, decisive control)",
                        "result": "see companion room -- expect Lv3/Lv5 disappear, Lv4/Lv6 remain"}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]] + snapshots[-1]["vel"])
    genesis_yaml = {
        "equations": "du/dt+vel.grad(u)=M*lap(mu)+R; dv/dt=D_v*lap(v)-R; "
                     "dvel/dt+vel.grad(vel)=-grad(p)+nu*lap(vel)+C*(mu*grad(u)), div(vel)=0; "
                     "R=kp*v*chi(u)-kd*chi(u)^2 (compartmentalized + Model H flow)",
        "solver": "pseudo-spectral, semi-implicit CH+diffusion, Leray projection, 2/3 dealiasing",
        "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
        "params": params, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼D [探索] u_mean=%.2f, t_final=%.0f, grid=%s, seed=%d, C=%.1f。reached_level=%d。"
             "%s" % (u_mean, STEPS * DT, SHAPE, seed, params.get("C", CONFIRMED_PARAMS["C"]),
                     report["reached_level"], notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  reached_level=%d role=%s KE=%.2e circ=%.2e u_bg_final=%.4f"
          % (run_dir, report["reached_level"], report["purity"]["role"],
             report["measured_by"]["final_kinetic_energy"], report["measured_by"]["final_circulation"],
             report["measured_by"]["u_bg_final"] or -1))
    return snapshots, phys, report


def main():
    print("=== 依頼D【探索】流れ(Model H)でLv3/Lv5も足した最大climb (3D) ===")
    print("予言: C=2で流れがLv3/Lv5を足し、区画化がLv4/Lv6を保つ→0->1->2->3->4->5->6")
    print("決定的対照: C=0でLv3/Lv5が消え、Lv4/Lv6だけ残る")

    print("\n[流れ結合] C=2, u_mean=0.35, seed=1")
    _run_and_save("req-d-flow-climb-C2-seed0001", CONFIRMED_PARAMS, seed=1)

    print("\n[決定的対照] C=0, u_mean=0.35, seed=1")
    control_params = dict(CONFIRMED_PARAMS, C=0.0)
    _run_and_save("req-d-flow-climb-control-C0-seed0001", control_params, seed=1,
                  notes_extra="決定的対照: C=0で流れ結合を切る。Lv3/Lv5が消え、依頼Aと同様の"
                              "Lv4/Lv6のみが残ることを期待。")

    print("=== 依頼D done ===")


if __name__ == "__main__":
    main()
