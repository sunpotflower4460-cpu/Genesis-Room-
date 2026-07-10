#!/usr/bin/env python3
"""genesis/g003_model_h_3d.py -- 依頼1: G003 Model H (3D)。Cahn-Hilliard 相分離場 phi と非圧縮
Navier-Stokes 流れ u が、毛管力 C*mu*grad(phi) だけで結合する（部品接続でなく一つの自由エネルギー
・化学ポテンシャル・応力場から）:

    d(phi)/dt + u.grad(phi) = M*lap(mu),  mu = phi^3 - phi - kappa*lap(phi)   (Cahn-Hilliard)
    d(u)/dt + u.grad(u) = -grad(p) + nu*lap(u) + C*(mu*grad(phi)),  div(u)=0  (Navier-Stokes)

完成した液滴・界面は一切初期条件に置かない：phi は一様平均 + 微小ノイズから開始する
(requests/...md 依頼1 §A initial_state)。

数値法：CH の biharmonic 剛性と粘性拡散は半陰的スペクトル（無条件安定）、移流・毛管力は
陽的擬スペクトル（2/3 則で dealiasing）。div(u)=0 は毎ステップ Leray 射影で強制する。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import emergence, io, diagnostics as diag  # noqa: E402
from genesis.solvers import k_grid, dealias_mask, leray_project  # noqa: E402

MODEL_ID = "g003_model_h_3d"

DEFAULTS = {"M": 1.0, "kappa": 1.0, "nu": 1.0, "C": 2.0}


def make_initial(shape, noise_amplitude, rng, mean=0.0):
    """一様平均 + 微小ノイズ。液滴/界面は入れない（NO seeded droplets/interfaces）。"""
    ndim = len(shape)
    phi = mean + noise_amplitude * rng.standard_normal(shape)
    u = [np.zeros(shape) for _ in range(ndim)]
    return phi, u


def _hat(field, dealias):
    return np.fft.fftn(field) * dealias


def step(phi, u, dt, p, kk, k2, dealias):
    """One semi-implicit pseudo-spectral step. Returns (phi_new, u_new, mu_real)."""
    ndim = phi.ndim
    phihat_full = np.fft.fftn(phi)  # undealiased, used for the linear/implicit part
    phihat = phihat_full * dealias  # dealiased copy for nonlinear products

    lap_phi = np.real(np.fft.ifftn(-k2 * phihat_full))
    mu = phi ** 3 - phi - p["kappa"] * lap_phi
    grad_phi = [np.real(np.fft.ifftn(1j * kk[a] * phihat)) for a in range(ndim)]

    advect_phi = sum(u[a] * grad_phi[a] for a in range(ndim))
    phi3_hat = _hat(phi ** 3, dealias)
    explicit_phi_hat = _hat(-advect_phi, dealias) - p["M"] * k2 * phi3_hat
    # implicit linear part: d(phihat)/dt ⊃ +M*k2*phihat - M*kappa*k2^2*phihat (see module docstring
    # derivation in requests/...md background); backward-Euler removes the k^4 stiffness.
    denom_phi = 1.0 - dt * p["M"] * k2 + dt * p["M"] * p["kappa"] * k2 ** 2
    phihat_new = (phihat_full + dt * explicit_phi_hat) / denom_phi
    phi_new = np.real(np.fft.ifftn(phihat_new))

    uhat = [_hat(ua, dealias) for ua in u]
    grad_u = [[np.real(np.fft.ifftn(1j * kk[b] * uhat[a])) for b in range(ndim)] for a in range(ndim)]
    advect_u = [sum(u[b] * grad_u[a][b] for b in range(ndim)) for a in range(ndim)]
    capillary = [p["C"] * mu * grad_phi[a] for a in range(ndim)]
    rhs_hat = [_hat(-advect_u[a] + capillary[a], dealias) for a in range(ndim)]
    uhat_full = [np.fft.fftn(ua) for ua in u]
    u_star_hat = [uhat_full[a] + dt * rhs_hat[a] for a in range(ndim)]
    u_visc_hat = [u_star_hat[a] / (1.0 + dt * p["nu"] * k2) for a in range(ndim)]
    u_proj_hat = leray_project(u_visc_hat, kk, k2)
    u_new = [np.real(np.fft.ifftn(uh)) for uh in u_proj_hat]

    return phi_new, u_new, mu


def free_energy(phi, kk, p):
    """F = integral[ 1/4 (phi^2-1)^2 + kappa/2 |grad phi|^2 ] -- should be non-increasing (dissipative
    gradient flow); mass mean(phi) is the exactly-conserved quantity (integrity's conservation_drift).
    """
    phihat = np.fft.fftn(phi)
    grad_phi = [np.real(np.fft.ifftn(1j * kk[a] * phihat)) for a in range(phi.ndim)]
    grad2 = sum(g ** 2 for g in grad_phi)
    return float(np.mean(0.25 * (phi ** 2 - 1) ** 2 + 0.5 * p["kappa"] * grad2))


def run(shape, steps, dt, seed, params=None, snapshot_every=None, noise_amplitude=0.01, mean=0.0):
    """t=0（一様+微小ノイズ、完成形なし）から中断なく発展させ、snapshots を返す。"""
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    rng = np.random.default_rng(seed)
    phi, u = make_initial(shape, noise_amplitude, rng, mean=mean)
    kk, k2 = k_grid(shape)
    dealias = dealias_mask(shape)
    snapshot_every = snapshot_every or max(1, steps // 12)

    f0 = free_energy(phi, kk, p)
    mass0 = float(np.mean(phi))
    snapshots = []
    diverged = False
    ke_series = []
    for t in range(steps):
        phi, u, mu = step(phi, u, dt, p, kk, k2, dealias)
        if not np.all(np.isfinite(phi)) or any(not np.all(np.isfinite(c)) for c in u):
            diverged = True
            break
        ke_series.append(float(0.5 * np.mean(sum(c ** 2 for c in u))))
        if t % snapshot_every == 0 or t == steps - 1:
            snapshots.append({"step": t, "field": phi.copy(), "u": [c.copy() for c in u]})
    f1 = free_energy(phi, kk, p)
    mass1 = float(np.mean(phi))
    phys = {
        "free_energy_drift": float(f1 - f0),
        "free_energy_monotonic_decrease": bool(f1 <= f0 + 1e-6),
        "mass_drift": float(mass1 - mass0),
        "diverged": diverged,
        "final_kinetic_energy": ke_series[-1] if ke_series else 0.0,
        "kinetic_energy_series": ke_series,
    }
    if not snapshots:
        snapshots = [{"step": 0, "field": phi.copy(), "u": [c.copy() for c in u]}]
    return snapshots, phys


# --- Room 依頼1 (G003) 実験計画: 検証 + 探索 + 決定的対照 + 解像度収束 ------------------

DT = 0.05
STEPS = 2000  # t_final = 100
SNAPSHOT_EVERY = 100
N_HI = 48
N_LO = 32
PRIMARY_PARAMS = dict(DEFAULTS)


def _coarsening_series(snapshots):
    times = [s["step"] * DT for s in snapshots]
    lengths = [diag.coarsening_length(s["field"]) for s in snapshots]
    return times, lengths


def _fit_power_law(times, lengths):
    """L(t) ~ t^n: least-squares fit of log(L) vs log(t), t>0 only. Returns n (None if <2 points)."""
    pts = [(t, l) for t, l in zip(times, lengths) if t > 0 and l > 0]
    if len(pts) < 2:
        return None
    logt = np.log([p[0] for p in pts])
    logl = np.log([p[1] for p in pts])
    return float(np.polyfit(logt, logl, 1)[0])


def _level_report_model_h(snapshots, phys, control_snapshots=None, control_phys=None):
    """compute_level_report は複素位相巻き(Lv2)前提なので、実数の phi 場向けに Lv2 を
    connected_components ベースへ差し替え、Lv5（co-differentiation）を C=0 対照との比較で判定する。
    """
    report = emergence.compute_level_report(snapshots, kind="model_h")
    final_phi = snapshots[-1]["field"]

    frac_pos = diag.occupied_fraction(final_phi > 0)
    n_pos, _, _ = diag.connected_components(final_phi > 0)
    n_neg, _, _ = diag.connected_components(final_phi <= 0)
    localization = bool(report["detected"]["difference"] and 0.0 < frac_pos < 1.0 and n_pos >= 1 and n_neg >= 1)
    report["detected"]["localization"] = localization
    report["measured_by"]["occupied_fraction_phase_pos"] = round(float(frac_pos), 6)
    report["measured_by"]["connected_components_pos"] = int(n_pos)
    report["measured_by"]["connected_components_neg"] = int(n_neg)

    reached = report["reached_level"]
    if localization:
        reached = max(reached, 2)
    spontaneous_motion = report["detected"].get("spontaneous_motion", False)
    if localization and spontaneous_motion:
        reached = max(reached, 3)

    co_occurring = bool(localization and spontaneous_motion)
    mutual_dependence = False
    if co_occurring and control_snapshots is not None and control_phys is not None:
        control_ke = control_phys.get("final_kinetic_energy", None)
        control_len = diag.coarsening_length(control_snapshots[-1]["field"])
        primary_len = diag.coarsening_length(final_phi)
        mutual_dependence = bool(control_ke is not None and control_ke < 1e-10 and primary_len > control_len)
        report["measured_by"]["control_C0_final_kinetic_energy"] = round(float(control_ke), 12)
        report["measured_by"]["control_C0_coarsening_length"] = round(float(control_len), 6)
        report["measured_by"]["primary_coarsening_length"] = round(float(primary_len), 6)
    report["detected"]["co_differentiation"] = bool(co_occurring and mutual_dependence)
    if report["detected"]["co_differentiation"]:
        reached = max(reached, 5)

    report["reached_level"] = reached
    report["candidate_level"] = min(reached + 1, 8)
    times, lengths = _coarsening_series(snapshots)
    report["measured_by"]["coarsening_rate_exponent"] = _fit_power_law(times, lengths)
    report["measured_by"]["coarsening_length_final"] = round(float(lengths[-1]), 6) if lengths else 0.0
    report["measured_by"]["final_kinetic_energy"] = round(float(phys["final_kinetic_energy"]), 10)
    role = "E" if reached >= 1 else "F"
    report["purity"] = {"per_object_labels": False, "external_optimum": False, "role": role}
    return report


def _run_and_save(room_id, shape, params, seed, mean=0.0, control_snapshots=None, control_phys=None,
                   notes_extra=""):
    snapshots, phys = run(shape, STEPS, DT, seed, params=params, snapshot_every=SNAPSHOT_EVERY, mean=mean)
    report = _level_report_model_h(snapshots, phys, control_snapshots, control_phys)

    resolutions_result = {"%dx%dx%d" % shape: report["measured_by"]["coarsening_length_final"]}
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"],
        resolutions_result=resolutions_result,
        seed_success={str(seed): report["reached_level"]},
        nan_or_clip=phys["diverged"],
    )
    integrity["free_energy_monotonic_decrease"] = phys["free_energy_monotonic_decrease"]
    integrity["free_energy_drift"] = phys["free_energy_drift"]

    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False,
        gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "C=0 (pure Cahn-Hilliard, no capillary coupling)",
                        "result": ("run separately; see NOTES for this room" if control_phys is None
                                   else "final_kinetic_energy=%.3e (u stays exactly 0 since forcing is "
                                        "identically 0 with C=0), coarsening_length=%.4f"
                                        % (control_phys["final_kinetic_energy"],
                                           diag.coarsening_length(control_snapshots[-1]["field"])))}],
    )

    checksum = io.checksum_of([snapshots[-1]["field"]] + snapshots[-1]["u"])
    genesis_yaml = {
        "equations": ("d(phi)/dt + u.grad(phi) = M*lap(mu), mu=phi^3-phi-kappa*lap(phi); "
                      "d(u)/dt + u.grad(u) = -grad(p) + nu*lap(u) + C*(mu*grad(phi)), div(u)=0"),
        "solver": "pseudo-spectral, semi-implicit CH/viscous, Leray projection, 2/3 dealiasing",
        "dt": DT, "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": params, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼1 [検証]: 2D で見つけた Lv1(分散/構造因子)->Lv2(界面/連結成分)->Lv3(KE>0)->"
             "Lv5(共分化、C=0対照で流れ消失+粗大化減速により相互依存を確認) が 3D 48^3 grid, "
             "t_final=%.0f (steps=%d, dt=%.3f) で立つかを検証。mode=coarse-global-3d 相当（依頼書の "
             "128^3 full-3d には未到達、下記注意点を参照）。%s"
             % (STEPS * DT, STEPS, DT, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output, figures={},
                                notes=notes)
    print("  wrote %s  reached_level=%d role=%s mass_drift=%.3e F_drift=%.4f"
          % (run_dir, report["reached_level"], report["purity"]["role"], phys["mass_drift"],
             phys["free_energy_drift"]))
    return snapshots, phys, report


def main():
    shape_hi = (N_HI, N_HI, N_HI)
    shape_lo = (N_LO, N_LO, N_LO)

    print("=== 依頼1 G003 Model H (3D) ===")
    print("[control] C=0 (pure CH, decisive control) seed=1 grid=%s" % (shape_hi,))
    control_params = dict(PRIMARY_PARAMS, C=0.0)
    control_snaps, control_phys, control_report = _run_and_save(
        "g003-model-h-3d-control-C0-seed0001", shape_hi, control_params, seed=1)

    print("[primary] default params, seeds=1,2,3, grid=%s" % (shape_hi,))
    primary_runs = {}
    for seed in (1, 2, 3):
        snaps, phys, report = _run_and_save(
            "g003-model-h-3d-seed%04d" % seed, shape_hi, PRIMARY_PARAMS, seed=seed,
            control_snapshots=control_snaps, control_phys=control_phys)
        primary_runs[seed] = (snaps, phys, report)

    print("[resolution] grid=%s vs %s, seed=1" % (shape_lo, shape_hi))
    _run_and_save("g003-model-h-3d-res32-seed0001", shape_lo, PRIMARY_PARAMS, seed=1,
                  control_snapshots=control_snaps, control_phys=control_phys,
                  notes_extra="解像度収束チェック（32^3 vs 48^3, 同一 seed/params/steps）。")

    print("[explore] off-critical mean=0.3 (droplet regime), seed=1")
    _run_and_save("g003-model-h-3d-offcritical-seed0001", shape_hi, PRIMARY_PARAMS, seed=1, mean=0.3,
                  control_snapshots=control_snaps, control_phys=control_phys,
                  notes_extra="[探索] off-critical (mean=0.3): critical(mean=0)のbicontinuousに対し"
                              "droplet形態になるかを比較。")

    print("[explore] high C=5, low nu=0.3 (push hydrodynamic instability)")
    push_params = dict(PRIMARY_PARAMS, C=5.0, nu=0.3)
    _run_and_save("g003-model-h-3d-highC-lownu-seed0001", shape_hi, push_params, seed=1,
                  control_snapshots=control_snaps, control_phys=control_phys,
                  notes_extra="[探索] 高C・低nu: 毛管力を強く/粘性を弱くして流体力学的不安定・"
                              "より速い粗大化が創発するか。")

    print("=== 依頼1 G003 done ===")


if __name__ == "__main__":
    main()
