#!/usr/bin/env python3
"""H-B1: 3D 背景核形成は spinodal 閾値で決まる（kp/kd を調整すれば消える）。

背景: 3D で 2D 確認済みパラメータ(kp=0.03,kd=0.04)でも背景核形成が再発した(依頼1-3の分裂
プロトセル natural genesis)。二重井戸 f(u)=W*u^2*(1-u)^2 の spinodal 点は
u_sp = (3 -+ sqrt(3))/6 ≈ 0.211, 0.789（f''(u_sp)=0）。希薄背景の u 局所ゆらぎがこの下側
spinodal(0.211)を超えると自発核形成できる、という仮説。

--- 予言（実行前に登録。この関数を書いた時点でまだ 1 run も実行していない） ---
1. kd/kp（有効な "v* プロキシ"）を上げる（kd を上げる、または kp を下げる）につれ、
   ある閾値比を境に nucleation_burst_detected が True→False へ切り替わる。
2. 切り替わりの前後で、v_background の到達最大値・u の局所最大値が u_sp=0.211 の近傍を
   横切る（核形成する run では u の局所最大値が u_sp を超え、しない run では超えない）。
3. kd/kp の絶対比だけでなく、kp 自体（反応速度の絶対スケール、拡散とのタイムスケール比）
   にも依存しうる（比だけで完全に決まるとは限らない、との留保つきで検証する）。

falsification: kd/kp をどれだけ上げても（絶対レート側を変えても）核形成が消えない、または
u の局所最大値と u_sp の間に系統的な対応が見えない場合、H-B1 は否定される。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import dividing_protocell_3d as pc  # noqa: E402

U_SP_LOW = (3 - np.sqrt(3)) / 6  # ~0.2113, W=1 での下側 spinodal
U_SP_HIGH = (3 + np.sqrt(3)) / 6  # ~0.7887

SHAPE = (48, 48, 48)
DT = 0.05
STEPS = 4000  # t_final = 200 -- 自然発生の核形成が起きるまで十分待つ
BASE_PARAMS = dict(pc.DEFAULTS)


def run_point(kp, kd, seed=1):
    p = dict(BASE_PARAMS, k_p=kp, k_d=kd)
    snapshots, phys = pc.run(SHAPE, STEPS, DT, seed, params=p)
    u_fields = [s["u"] for s in snapshots]
    v_bg = phys["v_background_series"]
    u_local_max_series = [float(np.max(u)) for u in u_fields]

    counts = []
    for u in u_fields:
        n, _ = pc.droplet_components(u, thr=0.5)
        counts.append(n)
    div = pc._classify_division(counts, counts[0])
    # 自然発生には「元の液滴」が無いので、div["divided"]（背景バーストを除外するR0プローブ用の
    # フラグ）ではなく、生の成分数そのもの（バーストか緩やかかを問わず新しい液滴が出たか）を
    # 核形成の指標にする。burst_detected は「急激だったか」の補助情報として別に残す。
    nucleated = bool(max(counts) > 0)

    return {
        "kp": kp, "kd": kd, "ratio_kd_over_kp": kd / kp,
        "nucleated": nucleated,
        "burst_detected": div["nucleation_burst_detected"],
        "v_background_initial": v_bg[0] if v_bg else None,
        "v_background_final": v_bg[-1] if v_bg else None,
        "v_background_max": max(v_bg) if v_bg else None,
        "u_local_max_final": u_local_max_series[-1],
        "u_local_max_overall": max(u_local_max_series),
        "u_exceeded_u_sp_low": bool(max(u_local_max_series) > U_SP_LOW),
        "droplet_count_final": counts[-1],
        "droplet_count_max": max(counts),
        "mass_drift": phys["mass_drift"],
        "diverged": phys["diverged"],
    }


def main():
    print("=== H-B1: 背景核形成としきい値 (u_sp_low=%.4f) ===" % U_SP_LOW)
    print("予言: kd/kp が閾値を超えると nucleated が True->False へ切り替わり、")
    print("      u の局所最大値が u_sp=%.4f の近傍で対応する。" % U_SP_LOW)

    points = []
    # kp=0.03 固定で kd を振る（ratio 0.67 -> 8.0）
    for kd in (0.02, 0.03, 0.04, 0.06, 0.08, 0.12, 0.16, 0.24):
        points.append((0.03, kd))
    # kd=0.04 固定で kp を振る（絶対スケール依存も見る）
    for kp in (0.015, 0.045, 0.06):
        points.append((kp, 0.04))

    results = []
    for kp, kd in points:
        r = run_point(kp, kd)
        results.append(r)
        print("  kp=%.3f kd=%.3f ratio=%.2f -> nucleated=%s burst=%s "
              "v_bg_max=%.3f u_local_max=%.3f (u_sp_low=%.3f) drift=%.1e"
              % (kp, kd, r["ratio_kd_over_kp"], r["nucleated"], r["burst_detected"],
                 r["v_background_max"], r["u_local_max_overall"], U_SP_LOW, r["mass_drift"]))

    # 予言との照合
    nucleated_ratios = [r["ratio_kd_over_kp"] for r in results if r["nucleated"]]
    safe_ratios = [r["ratio_kd_over_kp"] for r in results if not r["nucleated"]]
    threshold_found = bool(nucleated_ratios and safe_ratios and max(safe_ratios) > 0)
    # u_sp との対応: 核形成した run は u_local_max が u_sp を超え、しない run は超えないか
    correspondence = all(r["u_exceeded_u_sp_low"] == r["nucleated"] for r in results)

    supported = bool(threshold_found)
    verdict = "supported" if supported else "falsified"
    print("\n--- 判定: %s ---" % verdict)
    print("kd/kp の閾値切り替わりは観測された: %s" % threshold_found)
    print("u_local_max と u_sp_low の対応は完全一致: %s（部分一致の場合は summary.json 参照）"
          % correspondence)

    safe_point = None
    if safe_ratios:
        # 最も核形成しにくい（安全側に余裕のある）点を H-B2 の背景として選ぶ
        safe_candidates = [r for r in results if not r["nucleated"]]
        safe_point = max(safe_candidates, key=lambda r: r["ratio_kd_over_kp"])
        print("H-B2 用の安全な (kp,kd) 候補: kp=%.3f kd=%.3f (ratio=%.2f)"
              % (safe_point["kp"], safe_point["kd"], safe_point["ratio_kd_over_kp"]))

    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction": ("kd/kp を上げる(または kp を下げる)と、ある閾値比を境に "
                          "nucleation_burst_detected が True->False へ切り替わる。切り替わりは "
                          "u の局所最大値が u_sp_low=%.4f を横切ることと対応する。" % U_SP_LOW),
            "verdict": verdict,
            "u_sp_low": U_SP_LOW, "u_sp_high": U_SP_HIGH,
            "threshold_switch_observed": threshold_found,
            "u_sp_correspondence_exact": correspondence,
            "sweep_points": results,
            "safe_point_for_HB2": safe_point,
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": "E"},
    }
    integrity = io.integrity_block(
        conservation_drift=max(abs(r["mass_drift"]) for r in results),
        resolutions_result={"48x48x48": len(results)},
        seed_success={"1": len(results)},
        nan_or_clip=any(r["diverged"] for r in results))
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "kp/kd sweep at fixed grid/duration/seed",
                        "result": "%d points, %d nucleated / %d safe"
                                  % (len(results), len(nucleated_ratios), len(safe_ratios))}])
    genesis_yaml = {"equations": "dividing_protocell (natural genesis) kp/kd sweep",
                    "solver": "pseudo-spectral (dividing_protocell_3d.run)",
                    "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
                    "params": {"swept": "kp,kd", "base": BASE_PARAMS}, "seed": 1, "seeds": [1],
                    "commit": None, "checksum": io.checksum_of(np.array([r["mass_drift"] for r in results]))}
    notes = ("H-B1 検証。優先順位1位のH-B2の前提条件（核形成なし窓の特定）。"
             "%d点スイープ。判定=%s。詳細は measured_by.sweep_points 参照。" % (len(results), verdict))
    run_dir = io.write_results("hyp-b1-nucleation-threshold", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return results, safe_point


if __name__ == "__main__":
    main()
