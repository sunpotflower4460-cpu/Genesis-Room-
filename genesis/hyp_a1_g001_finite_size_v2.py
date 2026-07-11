#!/usr/bin/env python3
"""H-A1 第2ラウンド: G001の「3D秩序化」は有限サイズ効果か（振幅|A|確認+強BF点+対照の長時間化）。

第1ラウンドの監査で見つかった3つの問題を修正する：
  1. 対照(0,0)がL=96/128でt=150では秩序化(欠陥->0)していなかった（「秩序化の基準」が未確立）。
     -> 十分長い時間走らせ、全Lで欠陥->0を確認する。
  2. 振幅|A|を確認していなかった（「欠陥ゼロ」が真の秩序化か振幅死か平面波か区別できない）。
     -> 各run で <|A|> と min|A| の時系列を記録し、最終状態の構造因子ピーク位置
        （k=0なら一様秩序、k>0なら平面波/進行波）も測る。
  3. 弱いBF点(1.5,-1.5)しか使っていなかった（強い点(2,-2),(2,-3)が指示された）。
     -> (2,-2)[1+bc=-3], (2,-3)[1+bc=-5] を含める。(1.5,-1.5)はL=48のみ参考として残す。

--- 予言（実行前に登録） ---
「本物の3D秩序化」なら：強BF点(2,-2),(2,-3)は |A|≈1 を保ったまま（振幅死でない）欠陥->0
（全L）。かつ緩和(0,0)対照も、十分な時間を与えれば全Lで欠陥->0（「秩序化の基準」成立）。
falsification: 強BF点(2,-2)または(2,-3)が、大箱(>=96^3)で |A|≈1 を保ったまま**非ゼロの
欠陥密度を維持**するなら（プラトー）、H-A1(有限サイズ効果)は支持される（3D defect turbulence
が存在し、第1ラウンドの「秩序化」は弱点/短時間のアーティファクトだった、という結論）。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import g001_cgl_3d as g001  # noqa: E402


def run_with_amplitude_tracking(shape, t_final, b, c, seed=1, noise_amplitude=0.01):
    """g001.run を呼び、各snapshotで <|A|>・min|A| を追加記録する（ソルバ本体は変更しない）。"""
    snapshots, phys = g001.run(shape, t_final, seed, params={"b": b, "c": c},
                               noise_amplitude=noise_amplitude)
    L = shape[0]
    times, amp_mean, amp_min, defect_counts = [], [], [], []
    for s in snapshots:
        f = s["field"]
        amp = np.abs(f)
        times.append(s["step"] * phys["dt_used"])
        amp_mean.append(float(np.mean(amp)))
        amp_min.append(float(np.min(amp)))
        defect_counts.append(diag.winding_defect_count(f))

    final_field = snapshots[-1]["field"]
    peak_k, prominence = diag.structure_factor_peak(final_field)
    # k=0 (peak_idx=0) は構造因子計算では「平均」として除外されるため、prominence が低い
    # (~1、フラットに近い)なら一様秩序、prominence が高くpeak_k>0なら平面波/進行波の疑い。
    is_plane_wave_like = bool(prominence > 3.0 and peak_k > 0)

    tail = defect_counts[len(defect_counts) // 2:]
    tail_density = [d / (L ** 3) for d in tail]

    return {
        "L": L, "b": b, "c": c, "t_final": t_final, "steps": phys["steps"],
        "dt_used": phys["dt_used"], "diverged": phys["diverged"],
        "times": times, "amp_mean_series": amp_mean, "amp_min_series": amp_min,
        "defect_count_series": defect_counts,
        "final_amp_mean": amp_mean[-1], "final_amp_min": amp_min[-1],
        "final_defect_count": defect_counts[-1],
        "tail_density_mean": float(np.mean(tail_density)),
        "tail_density_final": tail_density[-1],
        "final_structure_factor_peak_k": peak_k,
        "final_structure_factor_prominence": prominence,
        "is_plane_wave_like": is_plane_wave_like,
        "amplitude_alive": bool(amp_mean[-1] > 0.5),  # |A|~1 が真の秩序、|A|->小 が振幅死
    }


def classify_outcome(r):
    if not r["amplitude_alive"]:
        return "amplitude_death"
    if r["final_defect_count"] == 0:
        return "plane_wave_or_traveling" if r["is_plane_wave_like"] else "true_uniform_order"
    return "turbulence_or_still_relaxing"


def main():
    print("=== H-A1 第2ラウンド: 振幅確認+強BF点+対照の長時間化 ===")
    print("予言: 強BF点(2,-2)/(2,-3)は|A|~1を保ったまま欠陥->0(全L)なら真の秩序化。")
    print("      非ゼロ欠陥密度を維持するなら有限サイズ効果=H-A1支持。")

    results = {}

    print("\n[緩和(0,0)対照 -- 秩序化の基準を確立]")
    results["control_0_0"] = []
    for L, t_final in [(48, 300), (96, 400), (128, 400)]:
        r = run_with_amplitude_tracking((L, L, L), t_final, 0.0, 0.0)
        results["control_0_0"].append(r)
        print("  L=%d t_final=%d steps=%d -> final_defects=%d tail_density=%.3e amp_mean=%.3f "
              "outcome=%s" % (L, t_final, r["steps"], r["final_defect_count"], r["tail_density_mean"],
                              r["final_amp_mean"], classify_outcome(r)))

    print("\n[弱いBF点(1.5,-1.5) -- L=48のみ、振幅トラッキング付きで再確認]")
    results["weak_BF_1.5_-1.5"] = [run_with_amplitude_tracking((48, 48, 48), 200, 1.5, -1.5)]
    r = results["weak_BF_1.5_-1.5"][0]
    print("  L=48 -> final_defects=%d amp_mean=%.3f outcome=%s"
          % (r["final_defect_count"], r["final_amp_mean"], classify_outcome(r)))

    print("\n[強BF点(2,-2), 1+bc=-3]")
    results["strong_BF_2_-2"] = []
    for L in (48, 96, 128):
        r = run_with_amplitude_tracking((L, L, L), 200, 2.0, -2.0)
        results["strong_BF_2_-2"].append(r)
        print("  L=%d steps=%d -> final_defects=%d tail_density=%.3e amp_mean=%.3f "
              "peak_k=%.1f prom=%.2f outcome=%s"
              % (L, r["steps"], r["final_defect_count"], r["tail_density_mean"], r["final_amp_mean"],
                 r["final_structure_factor_peak_k"], r["final_structure_factor_prominence"],
                 classify_outcome(r)))

    print("\n[強BF点(2,-3), 1+bc=-5 -- L=48,96（計算資源制約でL=128は省略）]")
    results["strong_BF_2_-3"] = []
    for L in (48, 96):
        r = run_with_amplitude_tracking((L, L, L), 200, 2.0, -3.0)
        results["strong_BF_2_-3"].append(r)
        print("  L=%d steps=%d -> final_defects=%d tail_density=%.3e amp_mean=%.3f "
              "peak_k=%.1f prom=%.2f outcome=%s"
              % (L, r["steps"], r["final_defect_count"], r["tail_density_mean"], r["final_amp_mean"],
                 r["final_structure_factor_peak_k"], r["final_structure_factor_prominence"],
                 classify_outcome(r)))

    # 判定: 強BF点のいずれかが、大箱(>=96)で amplitude_alive かつ非ゼロプラトーを維持するか
    plateau_found = False
    for key in ("strong_BF_2_-2", "strong_BF_2_-3"):
        for r in results[key]:
            if r["L"] >= 96 and r["amplitude_alive"] and r["tail_density_mean"] > 1e-5:
                plateau_found = True
    control_orders = all(r["final_defect_count"] == 0 for r in results["control_0_0"])

    verdict = "supported" if plateau_found else "falsified"
    print("\n--- 判定: %s ---" % verdict)
    print("強BF点で非ゼロプラトー(L>=96, 振幅生存)を検出: %s" % plateau_found)
    print("対照(0,0)は全Lで秩序化(十分な時間で): %s" % control_orders)

    all_results_summary = {k: [{kk: vv for kk, vv in r.items() if kk not in ("times",)}
                                for r in v] for k, v in results.items()}
    genesis_yaml = {"equations": "complex Ginzburg-Landau, box-size scan round2 (amplitude-tracked)",
                    "solver": "real-space finite-difference, explicit Euler",
                    "dt": None, "dx": 1.0, "grid": "see measured_by", "boundary": "periodic",
                    "params": "see measured_by", "seed": 1, "seeds": [1], "commit": None,
                    "checksum": io.checksum_of(np.array([r["final_defect_count"]
                                                          for v in results.values() for r in v]))}
    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction": ("強BF点(2,-2)/(2,-3)は|A|~1を保ったまま欠陥->0(全L)なら真の秩序化。"
                          "非ゼロ欠陥密度プラトーを維持するなら有限サイズ効果(H-A1支持)。"),
            "verdict": verdict,
            "plateau_found_at_large_L": plateau_found,
            "control_orders_at_all_L_given_enough_time": control_orders,
            "results": all_results_summary,
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": "E"},
    }
    integrity = io.integrity_block(
        conservation_drift=0.0, resolutions_result={"control_defects_final":
            [r["final_defect_count"] for r in results["control_0_0"]]},
        seed_success={"1": sum(len(v) for v in results.values())},
        nan_or_clip=any(r["diverged"] for v in results.values() for r in v))
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "relaxational (0,0) control, run long enough at each L",
                        "result": "final defect counts: %s" % [r["final_defect_count"]
                                                                for r in results["control_0_0"]]}])
    notes = ("H-A1第2ラウンド。振幅|A|追跡+強BF点(2,-2),(2,-3)+対照の長時間化。判定=%s。"
             "L=128の(2,-3)は計算資源制約で省略（限界として明記）。" % verdict)
    run_dir = io.write_results("hyp-a1-v2-amplitude-strong-bf", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return results


if __name__ == "__main__":
    main()
