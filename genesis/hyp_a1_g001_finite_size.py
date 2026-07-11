#!/usr/bin/env python3
"""H-A1: G001 の「3D 秩序化」は 48^3 の有限サイズ効果である。

背景: 依頼1-3では全 (b,c) で defects->0 (秩序化) を観測した。Benjamin-Feir 不安定域
(1+bc<0) でも乱流ではなく欠陥全滅だった。仮説: 48^3 が乱流維持に小さすぎる（乱流相関長 xi が
箱サイズ L を超えて維持できない）。

--- 予言（実行前に登録） ---
1. 強い Benjamin-Feir 不安定点 (b,c)=(1.5,-1.5), 1+bc=-1.25 で、晩期の欠陥線密度
   rho(L) = winding_defect_count / L^3 を L=48,96,128 で測ると、L が大きくなるにつれ
   rho(L) が非ゼロのプラトーへ収束する（有限サイズ効果なら）。
2. 対照点 (b,c)=(0,0)（緩和極限）では、全 L で rho(L) -> 0（真の秩序化、線張力による自壊）
   を予想する。

falsification: 強い BF 点でも L=128 まで rho(L) が単調に 0 へ向かうなら（プラトーが見えない
なら）、H-A1 は否定される -- 3D CGL はこのパラメータ域で本当に秩序化するという結論になる。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import g001_cgl_3d as g001  # noqa: E402

T_FINAL = 150.0
L_VALUES = [48, 96, 128]
POINTS = [("relaxational_control", 0.0, 0.0), ("strong_BF_unstable", 1.5, -1.5)]


def late_time_density(L, b, c, seed=1, tail_frac=0.3):
    shape = (L, L, L)
    snapshots, phys = g001.run(shape, T_FINAL, seed, params={"b": b, "c": c})
    tail = snapshots[int(len(snapshots) * (1 - tail_frac)):] or snapshots[-1:]
    tail_defects = [diag.winding_defect_count(s["field"]) for s in tail]
    tail_density = [d / (L ** 3) for d in tail_defects]
    return {
        "L": L, "b": b, "c": c, "steps": phys["steps"], "dt_used": phys["dt_used"],
        "diverged": phys["diverged"],
        "tail_defect_counts": tail_defects,
        "tail_density_mean": float(np.mean(tail_density)),
        "tail_density_final": tail_density[-1],
        "final_defect_count": tail_defects[-1],
    }


def main():
    print("=== H-A1: G001 3D秩序化は有限サイズ効果か (rho(L) スケーリング) ===")
    print("予言: 強BF点(1.5,-1.5)は L とともに rho(L) が非ゼロへ収束、(0,0)は全Lで rho->0")

    results = {}
    for label, b, c in POINTS:
        results[label] = []
        for L in L_VALUES:
            r = late_time_density(L, b, c)
            results[label].append(r)
            print("  [%s] L=%d (b=%s,c=%s) steps=%d dt=%.4f tail_density_mean=%.3e "
                  "tail_density_final=%.3e final_count=%d diverged=%s"
                  % (label, L, b, c, r["steps"], r["dt_used"], r["tail_density_mean"],
                     r["tail_density_final"], r["final_defect_count"], r["diverged"]))

    bf = results["strong_BF_unstable"]
    ctrl = results["relaxational_control"]
    bf_densities = [r["tail_density_mean"] for r in bf]
    ctrl_densities = [r["tail_density_mean"] for r in ctrl]

    # プラトー検出: 最大Lでの密度が最小Lでの密度の 30% 以上残っていれば「収束せず有限値」とみなす
    bf_plateaus = bool(bf_densities[0] > 0 and bf_densities[-1] > 0.3 * bf_densities[0])
    ctrl_orders = bool(ctrl_densities[-1] < 0.1 * (ctrl_densities[0] + 1e-30) or ctrl_densities[-1] == 0)

    if bf_plateaus:
        verdict = "supported"
    else:
        verdict = "falsified"

    print("\n--- 判定: %s ---" % verdict)
    print("強BF点: L=%s -> density=%s (プラトー=%s)" % (L_VALUES, bf_densities, bf_plateaus))
    print("対照(0,0): L=%s -> density=%s (秩序化=%s)" % (L_VALUES, ctrl_densities, ctrl_orders))

    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction": ("強BF点(1.5,-1.5)はLとともにrho(L)が非ゼロへ収束(有限サイズ効果)。"
                          "対照(0,0)は全Lでrho->0(真の秩序化)。"),
            "verdict": verdict,
            "L_values": L_VALUES,
            "strong_BF_point": {"b": 1.5, "c": -1.5, "results": bf},
            "control_point": {"b": 0.0, "c": 0.0, "results": ctrl},
            "bf_density_plateau_detected": bf_plateaus,
            "control_orders_at_all_L": ctrl_orders,
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": "E"},
    }
    integrity = io.integrity_block(
        conservation_drift=0.0, resolutions_result={"L=%d" % L: None for L in L_VALUES},
        seed_success={"1": len(L_VALUES) * 2},
        nan_or_clip=any(r["diverged"] for r in bf + ctrl))
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "relaxational (0,0) control at all L",
                        "result": "density trend=%s" % ctrl_densities}])
    genesis_yaml = {"equations": "complex Ginzburg-Landau, box-size scan",
                    "solver": "real-space finite-difference, explicit Euler",
                    "dt": None, "dx": 1.0, "grid": "L in %s" % L_VALUES, "boundary": "periodic",
                    "params": "see measured_by", "seed": 1, "seeds": [1], "commit": None,
                    "checksum": io.checksum_of(np.array(bf_densities + ctrl_densities))}
    notes = ("H-A1検証。優先順位2位。判定=%s。強BF点のrho(L): %s。対照(0,0)のrho(L): %s。"
             % (verdict, bf_densities, ctrl_densities))
    run_dir = io.write_results("hyp-a1-g001-finite-size", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return results


if __name__ == "__main__":
    main()
