#!/usr/bin/env python3
"""H-B2【最重要】: 清潔な 1->2 分裂は「臨界サイズ R_c すぐ上」でのみ起きる。

背景: 依頼1-3では大きい液滴（R0=15、臨界の数倍）でしか試しておらず、多断片化しか見なかった。
臨界サイズ R_c のすぐ上の狭い窓でしか清潔な 1->2 ピンチは起きない、という仮説。核形成の交絡は
H-B1 で特定した安全な (kp,kd) window を使って排除する。

--- 予言（実行前に登録） ---
安全な (kp,kd)（H-B1 の結果：kd/kp >= 2.0 で核形成なし。ここでは kp=0.03, kd=0.06 を使う）で、
単一 seed 液滴の R0 を細かく振ると：
  1. R0 < R_c        -> 安定（1 のまま、分裂なし、または縮んで消滅）
  2. R_c < R0 <~ 1.5*R_c -> 清潔な 1->2 分裂（娘サイズ比 -> 0.5 近傍、有意断片数=2、非対称でない）
  3. R0 >> R_c        -> 多断片化（依頼1-3で見たもの）

falsification: R_c 直上でも娘サイズ比が 0.5 に近づかず、清潔な二分裂がどの R0 でも出ないなら
H-B2 は否定される（分裂は本質的に非対称/多断片という結論）。

Phase 1（このファイル）: 粗いR0スイープで R_c を挟み込む。
Phase 2（hyp_b2_phase2.py、Phase1の結果を見てから書く）: R_c 直上を密に走らせる。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import dividing_protocell_3d as pc  # noqa: E402

SHAPE = (64, 64, 64)
DT = 0.05
STEPS = 5000  # t_final = 250 -- 核形成の交絡がないので長時間安全に走らせられる
SAFE_PARAMS = dict(pc.DEFAULTS, k_p=0.03, k_d=0.06)  # H-B1 で確認済みの安全点


def classify_fate(R0, seed=1, params=None, steps=STEPS):
    p = params or SAFE_PARAMS
    snapshots, phys = pc.run_droplet_probe(SHAPE, R0, steps, DT, seed, params=p)
    times, counts, size_lists, sig_counts = pc._analyze_components(snapshots, dx=1.0)
    div = pc._classify_division(sig_counts, sig_counts[0])
    initial_vol = size_lists[0][0] if size_lists and size_lists[0] else 0.0
    final_vol = sum(v for v in size_lists[-1] if v >= initial_vol * 0.05) if size_lists[-1] else 0.0

    if sig_counts[-1] == 0 and initial_vol > 0 and final_vol < initial_vol * 0.05:
        fate = "dissolved"
    elif div["divided"]:
        fate = "divided"
    else:
        fate = "stable"

    daughter_ratio = None
    if fate == "divided" and div["division_snapshot_index"] is not None:
        sizes_at_split = [v for v in size_lists[div["division_snapshot_index"]] if v >= initial_vol * 0.05]
        if len(sizes_at_split) >= 2:
            ordered = sorted(sizes_at_split, reverse=True)
            daughter_ratio = round(ordered[1] / ordered[0], 4) if ordered[0] > 0 else 0.0

    return {
        "R0": R0, "fate": fate, "initial_volume": initial_vol, "final_volume": final_vol,
        "significant_final_count": sig_counts[-1], "significant_max_count": div["max_count"],
        "daughter_ratio": daughter_ratio, "division_time": (
            times[div["division_snapshot_index"]] if div["division_snapshot_index"] is not None else None),
        "raw_final_count": counts[-1], "mass_drift": phys["mass_drift"], "diverged": phys["diverged"],
        "time_series_sig_counts": list(zip([round(t, 1) for t in times], sig_counts)),
    }


def main():
    print("=== H-B2 Phase 1: 粗い R0 スイープで R_c を挟み込む (kp=%.3f kd=%.3f, 安全window) ==="
          % (SAFE_PARAMS["k_p"], SAFE_PARAMS["k_d"]))
    print("予言: R0<R_c で安定、R_c<R0<~1.5Rcで清潔な1->2分裂、R0>>Rcで多断片化")

    R0_values = [4, 6, 8, 10, 12, 14, 16, 18]
    results = []
    for R0 in R0_values:
        r = classify_fate(R0)
        results.append(r)
        print("  R0=%2d -> fate=%-9s sig_final=%d sig_max=%d daughter_ratio=%s div_time=%s drift=%.1e"
              % (R0, r["fate"], r["significant_final_count"], r["significant_max_count"],
                 r["daughter_ratio"], r["division_time"], r["mass_drift"]))

    stable_R0 = [r["R0"] for r in results if r["fate"] in ("stable", "dissolved")]
    divided_R0 = [r["R0"] for r in results if r["fate"] == "divided"]
    R_c_lower = max(stable_R0) if stable_R0 else None
    R_c_upper = min(divided_R0) if divided_R0 else None

    print("\n粗いブラケット: 安定/消滅=%s, 分裂=%s" % (stable_R0, divided_R0))
    print("R_c はおよそ (%s, %s) の間" % (R_c_lower, R_c_upper))

    genesis_yaml = {"equations": "dividing_protocell droplet probe, R0 coarse sweep (safe window)",
                    "solver": "pseudo-spectral (dividing_protocell_3d.run_droplet_probe)",
                    "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
                    "params": SAFE_PARAMS, "seed": 1, "seeds": [1], "commit": None,
                    "checksum": io.checksum_of(np.array([r["mass_drift"] for r in results]))}
    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": False, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction": "R0<R_c: 安定。R_c<R0<~1.5Rc: 清潔な1->2分裂(比->0.5)。R0>>Rc: 多断片化。",
            "phase": "1 (coarse bracket)",
            "R0_sweep": results, "R_c_bracket": [R_c_lower, R_c_upper],
        },
        "purity": {"per_object_labels": True, "external_optimum": False, "role": "S"},
    }
    integrity = io.integrity_block(
        conservation_drift=max(abs(r["mass_drift"]) for r in results),
        resolutions_result={"64x64x64": len(results)}, seed_success={"1": len(results)},
        nan_or_clip=any(r["diverged"] for r in results))
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # R0 は調べたい量そのもの: 明示的に申告
        gate_encodes_conclusion_causality=False, gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "H-B1 safe (kp,kd) window -- no background nucleation confound",
                        "result": "kp=%.3f kd=%.3f (ratio=%.2f), confirmed nucleation-free in H-B1"
                                  % (SAFE_PARAMS["k_p"], SAFE_PARAMS["k_d"],
                                     SAFE_PARAMS["k_d"] / SAFE_PARAMS["k_p"])}])
    notes = ("H-B2 Phase1（粗いブラケット）。R_c は (%s, %s) の間と推定。Phase2で密に走らせる。"
             % (R_c_lower, R_c_upper))
    run_dir = io.write_results("hyp-b2-phase1-coarse-bracket", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return results, R_c_lower, R_c_upper


if __name__ == "__main__":
    main()
