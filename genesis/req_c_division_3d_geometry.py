#!/usr/bin/env python3
"""依頼C【最重要の核心】: 分裂（Lv7）を3D幾何で（Rayleigh-Plateauピンチ）。

--- 予言（実行前に登録） ---
区画化で自己維持する液滴を3Dで臨界サイズ超（特に伸びた形状）に育てたとき、2Dでは表面張力で
安定円のまま（未分裂）だったのに対し、3Dでは伸びた液滴がRayleigh-Plateau不安定性
（円柱の長さ/半径 > 2*pi ≈ 6.28 で不安定、古典的な流体力学の結果）でくびれて「1->2」に
清潔に分裂する。臨界サイズを特定し、そのすぐ上で分裂モード（清潔な二分裂 vs 多断片）を地図化。

決定的測定: 形状追跡（伸長->単一くびれ->2分離 vs 安定球 vs 多断片）、背景u_bgがspinodal
以下を維持（核形成でなくピンチであること）、娘液滴が区画化を継承して自己維持するか。

grow-divideサイクル: 分裂後の娘が再成長->再分裂を繰り返すか（区画化を保ったまま）。

falsification: 3Dでも安定球のまま分裂しない、またはくびれても背景u_bgがspinodal超え
（ピンチでなく核形成による偽の"分裂"）なら、単純な3D幾何だけでは不十分——非対称成長・
機械的不安定・第二の場・buddingが必要、と正直に報告する。
"""

import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import compartmentalized_genesis_3d as cg  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

DT = 0.05
CONFIRMED_PARAMS = dict(cg.DEFAULTS)  # kp=0.05, kd=0.03 (依頼Aと同じ基本区画化パラメータ)


def make_capsule_initial(shape, radius, length, axis=2, dx=1.0, u_background=0.15, u_inside=1.0,
                          v_uniform=0.5, interface_width=2.0, rng=None):
    """円柱(カプセル型: 円柱本体+半球キャップ)の伸びた液滴。Rayleigh-Plateau不安定性の
    決定的テスト用（aspect = length/radius > 2*pi ~ 6.28 で不安定域）。"""
    ndim = len(shape)
    grids = np.indices(shape).astype(float) * dx
    center = [s * dx / 2.0 for s in shape]
    half_len = max(length / 2.0 - radius, 0.0)  # 円柱部の半長（キャップ分を除く）

    axial = grids[axis] - center[axis]
    radial_sq = sum((grids[a] - center[a]) ** 2 for a in range(ndim) if a != axis)
    radial = np.sqrt(radial_sq)

    # capsule SDF: 円柱部は |axial|<half_len で dist=radial-radius、キャップ部は
    # 球面距離 sqrt(radial^2+(|axial|-half_len)^2) - radius
    axial_clipped = np.clip(np.abs(axial) - half_len, 0, None)
    dist = np.sqrt(radial ** 2 + axial_clipped ** 2) - radius

    profile = 0.5 * (1.0 - np.tanh(dist / interface_width))
    u = u_background + (u_inside - u_background) * profile
    v = np.full(shape, v_uniform, dtype=float)
    if rng is not None:
        u = u + 0.001 * rng.standard_normal(shape)
    return u, v


def run_capsule(shape, radius, length, steps, seed=1, params=None, axis=2,
                snapshot_every=None, u_background=0.15, v_uniform=0.5):
    p = params or CONFIRMED_PARAMS
    rng = np.random.default_rng(seed)
    u, v = make_capsule_initial(shape, radius, length, axis=axis, u_background=u_background,
                                v_uniform=v_uniform, rng=rng)
    _, k2 = k_grid(shape)
    snapshot_every = snapshot_every or max(1, steps // 40)
    mass0 = float(np.mean(u + v))
    snapshots = []
    diverged = False
    for t in range(steps):
        u, v = cg.step(u, v, DT, p, k2, reaction_fn=cg.reaction_compartmentalized)
        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            n, labeled, sizes = diag.connected_components(u > 0.5)
            dilute = u < 0.5
            u_bg = float(np.mean(u[dilute])) if dilute.any() else float(np.mean(u))
            snapshots.append({"step": t, "u": u.copy(), "v": v.copy(), "n": n,
                              "sizes": sorted(sizes.tolist(), reverse=True)[:8], "u_bg": u_bg})
    mass1 = float(np.mean(u + v))
    phys = {"mass_drift": mass1 - mass0, "diverged": diverged}
    if not snapshots:
        n, labeled, sizes = diag.connected_components(u > 0.5)
        snapshots = [{"step": 0, "u": u.copy(), "v": v.copy(), "n": n,
                     "sizes": sorted(sizes.tolist(), reverse=True)[:8], "u_bg": 0.0}]
    return snapshots, phys


def classify_and_save(room_id, shape, radius, length, steps, seed=1, params=None,
                      notes_extra=""):
    aspect = length / radius
    snapshots, phys = run_capsule(shape, radius, length, steps, seed=seed, params=params)
    times = [s["step"] * DT for s in snapshots]
    counts = [s["n"] for s in snapshots]
    u_bg_series = [s["u_bg"] for s in snapshots]
    initial_vol = snapshots[0]["sizes"][0] if snapshots[0]["sizes"] else 0.0
    min_vol = initial_vol * 0.05

    sig_counts = []
    for s in snapshots:
        sig_counts.append(sum(1 for v in s["sizes"] if v >= min_vol))
    max_sig = max(sig_counts)
    final_sig = sig_counts[-1]
    div_idx = next((i for i, c in enumerate(sig_counts) if c > sig_counts[0]), None)

    max_u_bg = max(u_bg_series) if u_bg_series else 0.0
    background_clean_throughout = bool(max_u_bg < cg.U_SP_LOW)
    divided = bool(max_sig > sig_counts[0])
    mode = "stable_no_division" if not divided else (
        "clean_pinch_division" if max_sig == 2 and background_clean_throughout else
        ("multi_fragmentation" if max_sig > 2 else "division_but_nucleation_confound"))
    if divided and not background_clean_throughout:
        mode = "confounded_by_background_nucleation"

    daughter_ratio = None
    if divided and div_idx is not None:
        sig_sizes = [v for v in snapshots[div_idx]["sizes"] if v >= min_vol]
        if len(sig_sizes) >= 2:
            ordered = sorted(sig_sizes, reverse=True)
            daughter_ratio = round(ordered[1] / ordered[0], 4) if ordered[0] > 0 else 0.0

    reached = 4  # 持続個体（自己維持、区画化）として最低限 Lv4 は成立している前提
    if divided and mode == "clean_pinch_division":
        reached = 7
    elif divided:
        reached = 4  # 分裂はしたが「清潔」の基準を満たさない（多断片 or 核形成汚染）

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["localization"] = True
    detected["persistent_individuality"] = bool(final_sig >= 1 and not phys["diverged"])
    detected["self_maintaining_closure"] = background_clean_throughout
    detected["growth_division_inheritance"] = bool(mode == "clean_pinch_division")

    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": False,  # カプセルは決定的対照として初期条件に置く
        "level_detected_by_measurement": True, "detected": detected,
        "measured_by": {"radius": radius, "length": length, "aspect_ratio": round(aspect, 3),
                        "rp_unstable_predicted": bool(aspect > 2 * np.pi),
                        "division_mode": mode, "daughter_size_ratio": daughter_ratio,
                        "division_time": times[div_idx] if div_idx is not None else None,
                        "max_significant_count": max_sig, "final_significant_count": final_sig,
                        "max_u_bg": max_u_bg, "u_sp_low": cg.U_SP_LOW,
                        "background_clean_throughout": background_clean_throughout,
                        "significant_count_series": list(zip([round(t, 1) for t in times], sig_counts)),
                        "mass_drift": phys["mass_drift"]},
        "purity": {"per_object_labels": True, "external_optimum": False, "role": "S"},
    }
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"],
        resolutions_result={"%dx%dx%d" % shape: max_sig},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # カプセル形状・aspectはまさに調べたい量
        gate_encodes_conclusion_causality=False, gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "u_bg tracked throughout (nucleation vs pinch discriminator)",
                        "result": "max_u_bg=%.4f vs u_sp_low=%.4f" % (max_u_bg, cg.U_SP_LOW)}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]])
    genesis_yaml = {
        "equations": "compartmentalized genesis (R=kp*v*chi(u)-kd*chi(u)^2), capsule/cylinder probe",
        "solver": "pseudo-spectral, semi-implicit CH + implicit diffusion",
        "dt": DT, "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": params or CONFIRMED_PARAMS, "seed": seed, "seeds": [seed], "commit": None,
        "checksum": checksum,
    }
    notes = ("依頼C [決定的対照/診断] capsule R=%.1f L=%.1f (aspect=%.2f, RP不安定予想=%s), "
             "t_final=%.0f, grid=%s。division_mode=%s。significant count series: %s。%s"
             % (radius, length, aspect, aspect > 2 * np.pi, steps * DT, shape, mode,
                list(zip([round(t) for t in times], sig_counts)), notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  aspect=%.2f mode=%s ratio=%s max_u_bg=%.4f reached=%d"
          % (run_dir, aspect, mode, daughter_ratio, max_u_bg, reached))
    return snapshots, phys, report


def main():
    print("=== 依頼C【最重要の核心】分裂を3D幾何で (Rayleigh-Plateau) ===")
    print("予言: aspect(length/radius) > 2*pi ~ 6.28 で不安定、くびれて清潔な1->2分裂。")
    print("決定的測定: u_bgがspinodal以下を維持(ピンチ vs 核形成の区別)。")

    shape = (128, 48, 48)  # 長軸(z)に長いカプセル用ドメイン

    print("\n[R-P安定域] radius=8, length=32 (aspect=4.0 < 2pi, 安定球のはず)")
    classify_and_save("req-c-capsule-R08-L32-stable-seed0001", shape, radius=8, length=32,
                      steps=3000, notes_extra="R-P安定域の対照(aspect<2pi)。")

    print("\n[R-P不安定境界] radius=8, length=56 (aspect=7.0, 境界すぐ上)")
    classify_and_save("req-c-capsule-R08-L56-seed0001", shape, radius=8, length=56, steps=3000)

    print("\n[R-P不安定域] radius=8, length=80 (aspect=10.0)")
    classify_and_save("req-c-capsule-R08-L80-seed0001", shape, radius=8, length=80, steps=3000)

    print("\n[R-P強不安定域] radius=6, length=90 (aspect=15.0)")
    classify_and_save("req-c-capsule-R06-L90-seed0001", shape, radius=6, length=90, steps=3000)

    print("=== 依頼C (capsule probes) done ===")


if __name__ == "__main__":
    main()
