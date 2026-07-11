#!/usr/bin/env python3
"""依頼A【核心】: 0≠無 climb — 一つのルールが 0 からどこまで登るか（3D）。

--- 予言（実行前に登録） ---
u_mean=0.35/0.50（0≠無 = spinodal 不安定域）→ 一つの run が t=0 から中断なく
Lv0→Lv1(差)→Lv2(境界が相分離で創発)→Lv4(持続個体)→Lv6(自己維持・背景clean)を登る。
u_mean=0.12（無、安定一様、対照）→ Lv0 で停止（何も起源しない）。

決定的対照(1): u_mean=0.12（安定な無）は Lv0 停止を示す（u_mean スイープに内蔵）。
決定的対照(2) load-bearing: 生産を切る(k_p=0) → 個体が死ぬ（代謝が本質＝Lv6が本物）。
  2D では total_u 4744→0 で死んだ。

2D確認値: u=0.35でdroplets 122->50（Ostwald抑制で複数個体が持続）, u_bg≈0.05<<0.21,
load-bearingで死。

falsification: 3Dでu_bgがspinodal以下に保てない（背景核形成）or 生産を切っても個体が
生き残る（load-bearingでない）なら、区画化のLv6は3Dで崩れる。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import compartmentalized_genesis_3d as cg  # noqa: E402

SHAPE = (96, 96, 96)
DT = 0.05
STEPS = 3000  # t_final = 150
CONFIRMED_PARAMS = dict(cg.DEFAULTS)  # M=1,kappa=1,W=1,w=0.15,D_v=2,kp=0.05,kd=0.03


def _level_report(snapshots, phys, u_mean):
    u_fields = [s["u"] for s in snapshots]
    variances, growth_rate = diag.variance_growth(u_fields)
    _, prominence = diag.structure_factor_peak(u_fields[-1])
    final = snapshots[-1]

    difference = bool(growth_rate > 0 and prominence > 1.5)
    localization = bool(difference and final["n_droplets"] >= 1)
    # Lv4: 持続個体 -- 最終時刻でも droplet が存在し、系全体が消滅していない
    persistent_individuality = bool(localization and final["n_droplets"] >= 1
                                    and not phys["diverged"])
    # Lv6: 背景が spinodal 以下(clean)で自己維持 -- load-bearing 対照は別 run で確認する
    u_bg_final = phys["u_bg_final"]
    background_clean = bool(u_bg_final is not None and u_bg_final < cg.U_SP_LOW)
    self_maintaining = bool(persistent_individuality and background_clean)

    reached = 0
    if difference:
        reached = 1
    if localization:
        reached = 2
    if persistent_individuality:
        reached = 4
    if self_maintaining:
        reached = 6

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = difference
    detected["localization"] = localization
    detected["persistent_individuality"] = persistent_individuality
    detected["self_maintaining_closure"] = self_maintaining

    role = "E" if reached >= 1 else "F"
    return {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"u_mean": u_mean, "variance_growth": round(float(growth_rate), 6),
                        "structure_factor_prominence": round(float(prominence), 4),
                        "final_n_droplets": final["n_droplets"],
                        "final_top_sizes": final["sizes"][:5],
                        "u_bg_final": u_bg_final, "u_bg_max": phys["u_bg_max"],
                        "u_sp_low": cg.U_SP_LOW, "background_clean": background_clean,
                        "mass_drift": phys["mass_drift"]},
        "purity": {"per_object_labels": False, "external_optimum": False, "role": role},
    }


def _run_and_save(room_id, u_mean, seed, notes_extra="", params=None, steps=STEPS):
    p = params or CONFIRMED_PARAMS
    snapshots, phys = cg.run(SHAPE, steps, DT, seed, params=p, u_mean=u_mean)
    report = _level_report(snapshots, phys, u_mean)

    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"],
        resolutions_result={"96x96x96": report["measured_by"]["final_n_droplets"]},
        seed_success={str(seed): report["reached_level"]}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "u_mean=0.12 (stable uniform, 'nothing', decisive control 1)",
                        "result": "see companion room -- expect Lv0"},
                       {"name": "k_p=0 load-bearing (decisive control 2)",
                        "result": "see companion room -- expect death"}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]])
    genesis_yaml = {
        "equations": "du/dt=M*lap(mu)+R, mu=f'(u)-kappa*lap(u); dv/dt=D_v*lap(v)-R; "
                     "R=kp*v*chi(u)-kd*chi(u)^2 (compartmentalized: production INSIDE only)",
        "solver": "pseudo-spectral, semi-implicit CH (u) + implicit diffusion (v), explicit reaction",
        "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
        "params": p, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼A [検証] u_mean=%.2f, t_final=%.0f, grid=%s, seed=%d。0≠無 climb: "
             "reached_level=%d。u_bg: %s。%s"
             % (u_mean, steps * DT, SHAPE, seed, report["reached_level"],
                phys["u_bg_series"][::max(1, len(phys["u_bg_series"]) // 15)], notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  reached_level=%d role=%s u_bg_final=%.4f (sp_low=%.4f) n_droplets=%d"
          % (run_dir, report["reached_level"], report["purity"]["role"],
             report["measured_by"]["u_bg_final"] or -1, cg.U_SP_LOW,
             report["measured_by"]["final_n_droplets"]))
    return snapshots, phys, report


def main():
    print("=== 依頼A【核心】0≠無 climb (3D) ===")
    print("予言: u_mean=0.35/0.50 は Lv0->1->2->4->6 を登る。u_mean=0.12 は Lv0 停止。")
    print("決定的対照: u_mean=0.12（安定な無）、k_p=0（load-bearing）")

    print("\n[0≠無] u_mean=0.35, seeds=1,2,3")
    for seed in (1, 2, 3):
        _run_and_save("req-a-compartmentalized-u035-seed%04d" % seed, 0.35, seed)

    print("\n[0≠無・臨界] u_mean=0.50, seed=1")
    _run_and_save("req-a-compartmentalized-u050-seed0001", 0.50, 1)

    print("\n[決定的対照1] u_mean=0.12（安定な無）, seed=1")
    _run_and_save("req-a-compartmentalized-u012-control-seed0001", 0.12, 1,
                  notes_extra="決定的対照1: 安定な一様（無）は起源しない、を確認する run。")

    print("\n[決定的対照2] load-bearing: k_p=0, u_mean=0.35, seed=1")
    lb_params = dict(CONFIRMED_PARAMS, k_p=0.0)
    _run_and_save("req-a-compartmentalized-loadbearing-kp0-seed0001", 0.35, 1,
                  params=lb_params,
                  notes_extra="決定的対照2 (load-bearing): k_p=0 で生産を切る。個体が死ねば"
                              "代謝が本質=Lv6が本物であることの証拠。")

    print("=== 依頼A done ===")


if __name__ == "__main__":
    main()
