#!/usr/bin/env python3
"""依頼A【核心】: 能動収縮furrowで清潔な1->2分裂（3D）。

--- 予言（実行前に登録） ---
伸長した区画化液滴（自己維持する膜付き個体、aspect=2.5、RP不安定閾値2*pi=6.28未満で
passiveには安定）に、長軸に垂直な平面のガウス帯（赤道の"環"）に力を集中する能動収縮
furrowを当てると、清潔な「1->2」（正確に2個・娘サイズ比>0.4で対称・背景u_bg<spinodal）
になる。furrow OFF（f_str=0）では1個のまま（aspect=2.5はRP不安定域未満なので分裂しない）。

決定的対照(1): furrow OFF（f_str=0）-> 伸長液滴は分裂しない（1個 or 緩和）。
決定的対照(2): passive多断片化（前回のR0=8長尺カプセル、furrowなし長時間）-> 多断片化
（8個）に戻る -- 能動furrowが「多断片でなく清潔な2分」を出すことを示す。

2D確認値: furrow強度0.08〜0.20で2個・比1.00・u_bg=0.010。furrowなしは1個。

falsification: 3Dで能動furrowを当てても多断片化する or 2分にならないなら、furrow機構は
3Dで不十分（別の能動機構が要る）。

furrowは"機構"であって"答えの埋め込み"でない: 長軸・中心は液滴自身のボクセル分布のPCA
（慣性主軸と等価）から毎ステップ動的に読み取る。座標や向きをハードコードしない。sigma
（furrowの幅）も液滴自身の軸方向の広がり(90%ile)から算出する -- 「n=2にせよ」という命令
ではなく、幾何を読んで一点に力を集中する収縮環のアナロジー。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import compartmentalized_genesis_3d as cg  # noqa: E402
from genesis.req_c_division_3d_geometry import make_capsule_initial  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

DT = 0.05
CONFIRMED_PARAMS = dict(cg.DEFAULTS, k_d=0.02, f_str=0.14)  # 指示書のparameters


def compute_furrow_act(u, p, t, T, X, Y, Z, sigma_frac=0.35, min_voxels=50):
    """液滴自身のボクセル分布からPCAで長軸・中心を毎ステップ動的に読み取り、
    赤道面（長軸に垂直、中心通過）のガウス帯で能動収縮力actを計算する。
    幾何は液滴自身の状態から読む -- 軸・中心・sigmaいずれもハードコードしない。"""
    mask = u > 0.5
    n = int(mask.sum())
    if n < min_voxels:
        return np.zeros_like(u), None
    mx, my, mz = X[mask], Y[mask], Z[mask]
    centroid = np.array([mx.mean(), my.mean(), mz.mean()])
    pts = np.stack([mx - centroid[0], my - centroid[1], mz - centroid[2]], axis=1)
    cov = (pts.T @ pts) / n
    eigvals, eigvecs = np.linalg.eigh(cov)  # ascending
    axis = eigvecs[:, -1]  # 最大分散方向 = 長軸（慣性主軸の最小モーメント軸と等価）
    s_masked = pts @ axis
    half_extent = float(np.percentile(np.abs(s_masked), 90))
    sigma = max(sigma_frac * half_extent, 1e-6)
    sx = (X - centroid[0]) * axis[0] + (Y - centroid[1]) * axis[1] + (Z - centroid[2]) * axis[2]
    furrow = np.exp(-(sx / sigma) ** 2)
    prog = min(1.0, t / (0.6 * T))
    c = cg.chi(u, p["w"])
    act = -p["f_str"] * (0.3 + prog) * furrow * c
    return act, {"axis": axis.tolist(), "centroid": centroid.tolist(),
                 "half_extent": half_extent, "sigma": sigma, "prog": prog}


def step_active(u, v, dt, p, k2, t, T, X, Y, Z, furrow_on=True):
    """act は u 方程式のみに入る（v 方程式は R のみ: dv/dt = D_v*lap(v) - R、指示書の式通り）。"""
    if furrow_on:
        act, geom = compute_furrow_act(u, p, t, T, X, Y, Z)
    else:
        act, geom = np.zeros_like(u), None
    R_base = cg.reaction_compartmentalized(u, v, p)
    R_total = R_base + act
    fprime = 2.0 * p["W"] * u * (1.0 - u) * (1.0 - 2.0 * u)
    fprime_hat = np.fft.fftn(fprime)
    R_total_hat = np.fft.fftn(R_total)
    R_base_hat = np.fft.fftn(R_base)
    uhat = np.fft.fftn(u)
    vhat = np.fft.fftn(v)
    denom_u = 1.0 + dt * p["M"] * p["kappa"] * k2 ** 2
    uhat_new = (uhat + dt * (-p["M"] * k2 * fprime_hat + R_total_hat)) / denom_u
    denom_v = 1.0 + dt * p["D_v"] * k2
    vhat_new = (vhat - dt * R_base_hat) / denom_v
    return np.real(np.fft.ifftn(uhat_new)), np.real(np.fft.ifftn(vhat_new)), geom


def run_active(shape, radius, length, steps, dt, seed, params=None, axis=2, furrow_on=True,
               snapshot_every=None, u_background=0.15, v_uniform=0.5):
    p = dict(CONFIRMED_PARAMS)
    if params:
        p.update(params)
    rng = np.random.default_rng(seed)
    u, v = make_capsule_initial(shape, radius, length, axis=axis, u_background=u_background,
                                v_uniform=v_uniform, rng=rng)
    _, k2 = k_grid(shape)
    X, Y, Z = np.indices(shape).astype(float)
    T = steps * dt
    snapshot_every = snapshot_every or max(1, steps // 60)
    mass0 = float(np.mean(u + v))
    snapshots = []
    diverged = False
    geom_series = []
    for t in range(steps):
        u, v, geom = step_active(u, v, dt, p, k2, t * dt, T, X, Y, Z, furrow_on=furrow_on)
        if geom is not None:
            geom_series.append({"t": round(t * dt, 2), "prog": geom["prog"],
                                "half_extent": geom["half_extent"], "sigma": geom["sigma"]})
        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            n, _, sizes = diag.connected_components(u > 0.5)
            dilute = u < 0.5
            u_bg = float(np.mean(u[dilute])) if dilute.any() else float(np.mean(u))
            snapshots.append({"step": t, "u": u.copy(), "v": v.copy(), "n": n,
                              "sizes": sorted(sizes.tolist(), reverse=True)[:8], "u_bg": u_bg})
    mass1 = float(np.mean(u + v))
    if not snapshots:
        n, _, sizes = diag.connected_components(u > 0.5)
        snapshots = [{"step": 0, "u": u.copy(), "v": v.copy(), "n": n,
                     "sizes": sorted(sizes.tolist(), reverse=True)[:8], "u_bg": 0.0}]
    phys = {"mass_drift": mass1 - mass0, "diverged": diverged, "geom_series": geom_series}
    return snapshots, phys


def classify_and_save(room_id, shape, radius, length, steps, seed=1, params=None, furrow_on=True,
                      notes_extra=""):
    aspect = length / radius
    snapshots, phys = run_active(shape, radius, length, steps, DT, seed, params=params,
                                 furrow_on=furrow_on)
    times = [s["step"] * DT for s in snapshots]
    u_bg_series = [s["u_bg"] for s in snapshots]
    initial_vol = snapshots[0]["sizes"][0] if snapshots[0]["sizes"] else 0.0
    min_vol = initial_vol * 0.05

    sig_counts = [sum(1 for v in s["sizes"] if v >= min_vol) for s in snapshots]
    max_sig = max(sig_counts)
    final_sig = sig_counts[-1]
    div_idx = next((i for i, c in enumerate(sig_counts) if c > sig_counts[0]), None)

    max_u_bg = max(u_bg_series) if u_bg_series else 0.0
    background_clean_throughout = bool(max_u_bg < cg.U_SP_LOW)
    divided = bool(final_sig > sig_counts[0])
    mode = "stable_no_division" if not divided else (
        "clean_two_way_division" if final_sig == 2 and background_clean_throughout else
        ("multi_fragmentation" if final_sig > 2 else "division_but_nucleation_confound"))
    if divided and not background_clean_throughout:
        mode = "confounded_by_background_nucleation"

    daughter_ratio = None
    final_sizes = [v for v in snapshots[-1]["sizes"] if v >= min_vol]
    if final_sig >= 2:
        ordered = sorted(final_sizes, reverse=True)
        daughter_ratio = round(ordered[1] / ordered[0], 4) if ordered[0] > 0 else 0.0

    reached = 6  # 入力自体が自己維持する区画化液滴（Lv6相当）
    if divided and mode == "clean_two_way_division":
        reached = 7
    elif divided:
        reached = 6  # 分裂はしたが「清潔な2分」の基準を満たさない

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["localization"] = True
    detected["persistent_individuality"] = bool(final_sig >= 1 and not phys["diverged"])
    detected["self_maintaining_closure"] = background_clean_throughout
    detected["growth_division_inheritance"] = bool(mode == "clean_two_way_division")

    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": False,  # 伸長カプセルは決定的テスト用の初期条件
        "level_detected_by_measurement": True, "detected": detected,
        "measured_by": {"radius": radius, "length": length, "aspect_ratio": round(aspect, 3),
                        "furrow_on": furrow_on, "f_str": (params or CONFIRMED_PARAMS).get("f_str"),
                        "division_mode": mode, "daughter_size_ratio": daughter_ratio,
                        "division_time": times[div_idx] if div_idx is not None else None,
                        "max_significant_count": max_sig, "final_significant_count": final_sig,
                        "max_u_bg": max_u_bg, "u_sp_low": cg.U_SP_LOW,
                        "background_clean_throughout": background_clean_throughout,
                        "significant_count_series": list(zip([round(t, 1) for t in times], sig_counts)),
                        "mass_drift": phys["mass_drift"],
                        "geom_series_sampled": phys["geom_series"][::max(1, len(phys["geom_series"]) // 15)]
                        if phys["geom_series"] else []},
        "purity": {"per_object_labels": True, "external_optimum": False, "role": "S"},
    }
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"],
        resolutions_result={"%dx%dx%d" % shape: max_sig},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # カプセル形状・aspectは調べたい量そのもの
        gate_encodes_conclusion_causality=False,  # furrowはchi(u)ゲート+動的PCA軸読取のみ、n=2を知らない
        gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "furrow OFF (f_str=0, decisive control 1)",
                        "result": "see companion room -- expect no division"},
                       {"name": "u_bg tracked throughout (nucleation vs pinch discriminator)",
                        "result": "max_u_bg=%.4f vs u_sp_low=%.4f" % (max_u_bg, cg.U_SP_LOW)}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]])
    genesis_yaml = {
        "equations": "du/dt=M*lap(mu)+R+act, mu=f'(u)-kappa*lap(u); R=kp*v*chi(u)-kd*chi(u)^2; "
                     "act=-f_str*(0.3+prog)*exp(-(s/sigma)^2)*chi(u), s=coord along dynamic "
                     "(PCA-derived) long axis, prog=min(1,t/(0.6T))",
        "solver": "pseudo-spectral, semi-implicit CH + implicit diffusion, explicit furrow force",
        "dt": DT, "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": params or CONFIRMED_PARAMS, "seed": seed, "seeds": [seed], "commit": None,
        "checksum": checksum,
    }
    notes = ("依頼A [検証] active_division_3d capsule R=%.1f L=%.1f (aspect=%.2f), furrow_on=%s, "
             "t_final=%.0f, grid=%s。division_mode=%s, daughter_ratio=%s。significant count "
             "series: %s。%s"
             % (radius, length, aspect, furrow_on, steps * DT, shape, mode, daughter_ratio,
                list(zip([round(t) for t in times], sig_counts)), notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  furrow_on=%s mode=%s ratio=%s max_u_bg=%.4f reached=%d"
          % (run_dir, furrow_on, mode, daughter_ratio, max_u_bg, reached))
    return snapshots, phys, report


def main():
    print("=== 依頼A【核心】能動furrowで清潔な1->2分裂 (3D) ===")
    print("予言: aspect=2.5の伸長液滴、furrow ON->清潔な1->2。furrow OFF->1個のまま。")

    shape = (96, 96, 96)
    radius, length = 10, 25  # aspect=2.5 (RP不安定閾値2pi=6.28未満、passiveに安定なはず)
    steps = 3000  # t_final=150, prog深化は t=90 まで、その後60単位で分離を見る

    print("\n[本実験] furrow ON, seeds=1,2,3")
    for seed in (1, 2, 3):
        classify_and_save("req-a2-active-furrow-on-seed%04d" % seed, shape, radius, length,
                          steps, seed=seed, furrow_on=True)

    print("\n[決定的対照1] furrow OFF (f_str=0)")
    off_params = dict(CONFIRMED_PARAMS, f_str=0.0)
    classify_and_save("req-a2-active-furrow-off-control-seed0001", shape, radius, length,
                      steps, seed=1, params=off_params, furrow_on=False,
                      notes_extra="決定的対照1: furrow OFF (f_str=0)。分裂しない(1個)ことを期待。")

    print("=== 依頼A (active furrow) done ===")


if __name__ == "__main__":
    main()
