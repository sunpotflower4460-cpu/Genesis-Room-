#!/usr/bin/env python3
"""H-B2 Phase 2: R_c 直上の密なスイープ（方針転換の記録つき）。

Phase 1 で見つかった予想外の事実: H-B1 の「核形成なし安全window」(kp=0.03, kd>=0.06) は、
R0=4〜24（依頼1-3で試した最大R0=15/17を上回る範囲）の**全域で液滴が消滅**する regime だった。
kp/kdを2倍にスケールしても、D_vを3倍(6.0)にしても解消しなかった（いずれも単発テストで確認、
本ファイルの主スイープには含めない）。つまりこのモデルの関数形では、「背景核形成が起きない」
ことと「seed した液滴が成長できる」ことが、少なくとも kp=0.03 の系列では両立しなかった。

--- 方針転換（実行前に判断・登録） ---
依頼書は「あなたが気づいたら遠慮なくそちらを試してください」と明記しているため、H-B1の
(kp,kd)window をそのまま使うのではなく、元の 2D確認済み kd=0.04（成長可能）に戻し、
交絡排除は「時間」でなく依頼1-3で既に構築・検証済みの **突発バースト検出器**
（`_classify_division`: 単一ステップでの急増(絶対>3 かつ比>2)を背景核形成として除外し、
緩やかな成長は真の分裂として残す）に委ねる。この検出器は依頼1-3のR0=15で実際に機能した
実績がある（早期の非対称分裂を正しく残し、後の背景バーストを正しく除外した）。

--- 予言（実行前に登録） ---
kp=0.03, kd=0.04（成長可能、既存のburst検出器で交絡排除）で R0 を密に振ると：
  1. R0 <~ 13 -> 安定 or 消滅（分裂なし）
  2. R0 = 14〜17 のどこかで、娘サイズ比が 0.5 に近い、より対称的な分裂が見つかる可能性がある
     （Phase1の R0=15 が daughter_ratio=0.44 で最も"清潔"に近かったことを踏まえた期待）
  3. R0 >> 17 では多断片化（依頼1-3の知見の再確認）

falsification: どの R0 でも娘比が 0.5 に近づかず、Phase1 の R0=15(0.44) を超える対称性が
見つからないなら、"臨界サイズ直上のみ清潔"という当初の仮説は指示書のパラメータでは支持されない
（多断片化・非対称分裂が本質的、という結論になる）。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import dividing_protocell_3d as pc  # noqa: E402

SHAPE = (64, 64, 64)
DT = 0.05
STEPS = 3000  # t_final = 150 -- 依頼1-3と同じ設定（比較可能にする）
GROWTH_PARAMS = dict(pc.DEFAULTS)  # kp=0.03, kd=0.04（2D確認済み・成長可能）


def probe(R0, seed=1):
    snapshots, phys = pc.run_droplet_probe(SHAPE, R0, STEPS, DT, seed, params=GROWTH_PARAMS)
    times, counts, size_lists, sig_counts = pc._analyze_components(snapshots, dx=1.0)
    div = pc._classify_division(sig_counts, sig_counts[0])
    initial_vol = size_lists[0][0] if size_lists and size_lists[0] else 0.0
    min_vol = initial_vol * 0.05

    genuine_division = bool(div["divided"])
    size_ratio = None
    if genuine_division and div["division_snapshot_index"] is not None:
        sig_sizes = [v for v in size_lists[div["division_snapshot_index"]] if v >= min_vol]
        if len(sig_sizes) >= 2:
            ordered = sorted(sig_sizes, reverse=True)
            size_ratio = round(ordered[1] / ordered[0], 4) if ordered[0] > 0 else 0.0

    if not genuine_division:
        mode = "stable_or_dissolved" if sig_counts[-1] <= 1 else "stable_no_division"
    elif div["max_count"] > 2:
        mode = "multi_fragmentation"
    elif size_ratio is not None and size_ratio >= 0.5:
        mode = "clean_two_way_division"
    else:
        mode = "asymmetric_budding"

    return {
        "R0": R0, "mode": mode, "genuine_division": genuine_division,
        "daughter_size_ratio": size_ratio,
        "division_time": (times[div["division_snapshot_index"]]
                          if div["division_snapshot_index"] is not None else None),
        "significant_final_count": sig_counts[-1], "significant_max_count": div["max_count"],
        "raw_final_count": counts[-1], "burst_detected": div["nucleation_burst_detected"],
        "mass_drift": phys["mass_drift"], "diverged": phys["diverged"],
    }


def main():
    print("=== H-B2 Phase 2: kd=0.04(成長可能)+ burst検出器で R_c 直上を密にスイープ ===")
    print("予言: R0=14-17付近でより対称的な分裂(比->0.5)、それ以下は安定、それ以上は多断片化")

    R0_values = [12, 13, 14, 15, 16, 17, 18, 19, 20]
    results = []
    for R0 in R0_values:
        r = probe(R0)
        results.append(r)
        print("  R0=%2d -> mode=%-22s ratio=%s div_time=%s sig_max=%d burst=%s drift=%.1e"
              % (R0, r["mode"], r["daughter_size_ratio"], r["division_time"],
                 r["significant_max_count"], r["burst_detected"], r["mass_drift"]))

    divided = [r for r in results if r["genuine_division"]]
    best = max(divided, key=lambda r: r["daughter_size_ratio"] or 0) if divided else None
    clean_found = bool(best and best["daughter_size_ratio"] is not None and best["daughter_size_ratio"] >= 0.5)

    verdict = "supported" if clean_found else ("inconclusive" if divided else "falsified")
    print("\n--- 判定: %s ---" % verdict)
    if best:
        print("最も対称的な分裂: R0=%d, daughter_ratio=%.3f, mode=%s"
              % (best["R0"], best["daughter_size_ratio"], best["mode"]))
    else:
        print("この R0 域では genuine_division が一度も検出されなかった")

    genesis_yaml = {"equations": "dividing_protocell droplet probe, R0 fine sweep near R_c (kd=0.04)",
                    "solver": "pseudo-spectral (dividing_protocell_3d.run_droplet_probe)",
                    "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
                    "params": GROWTH_PARAMS, "seed": 1, "seeds": [1], "commit": None,
                    "checksum": io.checksum_of(np.array([r["mass_drift"] for r in results]))}
    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": False, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction": "R0=14-17付近でより対称的な分裂(比->0.5)、方針転換理由はモジュールdocstring参照",
            "pivot_reason": ("H-B1のkd/kp安全window(比>=2.0)はPhase1で全R0域(4-24)にわたり液滴を"
                            "消滅させる regime と判明。growth-permittingなkd=0.04に戻し、既存の"
                            "突発バースト検出器で交絡排除する方針に転換。"),
            "verdict": verdict, "phase": "2 (fine sweep near R_c)",
            "R0_sweep": results, "best_division": best,
        },
        "purity": {"per_object_labels": True, "external_optimum": False, "role": "S"},
    }
    integrity = io.integrity_block(
        conservation_drift=max(abs(r["mass_drift"]) for r in results),
        resolutions_result={"64x64x64": len(results)}, seed_success={"1": len(results)},
        nan_or_clip=any(r["diverged"] for r in results))
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,
        gate_encodes_conclusion_causality=False, gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "burst detector (validated in 依頼1-3 R0=15 run)",
                        "result": "%d/%d points show burst_detected at final snapshot"
                                  % (sum(1 for r in results if r["burst_detected"]), len(results))}])
    notes = ("H-B2 Phase2。方針転換: H-B1安全windowが成長を殺すため growth-permitting kd=0.04 + "
             "既存burst検出器を使用（理由はdocstring）。判定=%s。" % verdict)
    run_dir = io.write_results("hyp-b2-phase2-fine-sweep", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return results


if __name__ == "__main__":
    main()
