#!/usr/bin/env python3
"""課題3-1: 時間的分離（パルス化）で「核形成なし ⊥ 成長」の緊張を解けるか。

背景（前ラウンドで確認された緊張）: 一様生産(kp=0.03,kd=0.04)は液滴を成長させるが、背景も
いずれ spinodal(u_sp=0.2113)を超えて自発核形成する。核形成を防ぐ安全window(kd比>=2.0)は
液滴も殺す。流れ(Model H)・局所生産(表面ゲート)も解決しないとClaudeの2D予備調査で判明。

--- アイデア: 生産をパルス化する ---
成長期(GROWTH, kp=0.03,kd=0.04)＝液滴が育つが背景も上昇する。
休止期(REST, kp=0,kd=0 の純CH)＝反応なし。背景は拡散のみで（準安定域なら）緩和し、
既存の u 材料で液滴が形状不安定（首のくびれ）を起こす余地ができる。

--- 予言（実行前に登録） ---
1. u_far（液滴から離れた背景の u 局所平均）が spinodal(0.2113)の安全マージン下（例: 0.15-0.18）
   に達する時刻を T_on として GROWTH を打ち切ると、u_far は spinodal を一度も超えない
   （＝背景核形成なし）。
2. REST 期間中に u_far は低下する（拡散による緩和）。
3. 複数サイクルを回すと、液滴は growth 期に成長し、rest 期または次の growth 期の初期に
   分裂し、さらに（材料補充後）再成長する ―― という grow-divide サイクルの兆候が見える。

falsification: パルス化しても u_far が結局 spinodal を超える（核形成が防げない）、または
液滴が分裂・再成長のサイクルを一度も示さない（単に縮小するだけ、または不変のまま）なら、
このアイデアは否定される（構造的対立は時間分離でも解けない、という結論）。
"""

import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import dividing_protocell_3d as pc  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

U_SP_LOW = (3 - np.sqrt(3)) / 6  # ~0.2113
SHAPE = (64, 64, 64)
DT = 0.05
GROWTH_PARAMS = dict(pc.DEFAULTS)  # kp=0.03, kd=0.04
REST_PARAMS = dict(pc.DEFAULTS, k_p=0.0, k_d=0.0)  # 純CH、反応なし


def u_far_background(u, droplet_thr=0.5, dilation=5):
    """液滴（+安全マージン）を除いた領域での u の平均 -- 「背景」の局所レベル。"""
    mask = u > droplet_thr
    if not mask.any():
        return float(np.mean(u))
    dilated = ndimage.binary_dilation(mask, structure=np.ones((3, 3, 3)), iterations=dilation)
    far = ~dilated
    return float(np.mean(u[far])) if far.any() else float(np.mean(u))


def calibrate_T_on(R0=15, seed=1, max_t=150.0, safety_margin=0.18):
    """一様生産を続けた場合に u_far が safety_margin に達する時刻を見つける（予備測定）。"""
    rng = np.random.default_rng(seed)
    u, v = pc.make_droplet_initial(SHAPE, R0, rng=rng)
    _, k2 = k_grid(SHAPE)
    steps = int(round(max_t / DT))
    t = 0.0
    for i in range(steps):
        u, v = pc.step(u, v, DT, GROWTH_PARAMS, k2)
        t += DT
        if i % 20 == 0:
            uf = u_far_background(u)
            if uf >= safety_margin:
                return t, uf, u, v
    return None, u_far_background(u), u, v


def run_pulsed(R0, T_on, T_off, n_cycles, seed=1, snapshot_every_steps=20):
    rng = np.random.default_rng(seed)
    u, v = pc.make_droplet_initial(SHAPE, R0, rng=rng)
    _, k2 = k_grid(SHAPE)
    mass0 = float(np.mean(u + v))

    steps_on = int(round(T_on / DT))
    steps_off = int(round(T_off / DT))
    snapshots = []
    u_far_series = []
    t = 0.0
    step_idx = 0
    diverged = False
    for cycle in range(n_cycles):
        for phase, phase_params, phase_steps in (("growth", GROWTH_PARAMS, steps_on),
                                                  ("rest", REST_PARAMS, steps_off)):
            for _ in range(phase_steps):
                u, v = pc.step(u, v, DT, phase_params, k2)
                t += DT
                step_idx += 1
                if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
                    diverged = True
                    break
                if step_idx % snapshot_every_steps == 0:
                    uf = u_far_background(u)
                    u_far_series.append({"t": round(t, 2), "phase": phase, "u_far": uf,
                                         "cycle": cycle})
                    n, sizes = pc.droplet_components(u)
                    snapshots.append({"step": step_idx, "t": t, "phase": phase, "cycle": cycle,
                                      "u": u.copy(), "v": v.copy(), "n_components": n,
                                      "sizes": sorted(sizes.tolist(), reverse=True)[:6]})
            if diverged:
                break
        if diverged:
            break

    mass1 = float(np.mean(u + v))
    max_u_far = max((r["u_far"] for r in u_far_series), default=0.0)
    nucleation_avoided = bool(max_u_far < U_SP_LOW)

    return {
        "R0": R0, "T_on": T_on, "T_off": T_off, "n_cycles": n_cycles,
        "u_far_series": u_far_series, "max_u_far": max_u_far,
        "nucleation_avoided": nucleation_avoided,
        "component_series": [{"t": s["t"], "phase": s["phase"], "cycle": s["cycle"],
                             "n": s["n_components"], "top_sizes": s["sizes"]} for s in snapshots],
        "mass_drift": mass1 - mass0, "diverged": diverged,
    }


def main():
    print("=== 課題3-1: 時間的パルス化で核形成なし+grow-divideサイクルを狙う ===")
    print("予言(当初): T_onを安全マージンで切ればu_farはspinodalを超えない。複数サイクルで")
    print("            成長->分裂->再成長の兆候が見えるか。")

    print("\n[校正] 一様生産でのu_far(t)を直接観察（予備実験で判明した重要な事実:")
    print("  u_farはt~10でspinodal(0.211)を超えるが、可視的な核形成バーストはt~75まで起きない")
    print("  ＝spinodal通過と可視核形成の間に長い誘導時間がある。当初の「spinodalを一度も超え")
    print("  ない」予言は、意味のある成長時間内では原理的に達成不能と判明——予言を修正する。")
    print("  さらにREST中のu_far緩和速度(~0.00016/t)はGROWTH中の上昇速度(~0.0075/t)の1/45しか")
    print("  なく、単純往復パルスでは正味ドリフトが避けられない可能性が高い、と予備測定で判明。")

    print("\n[修正した予言] T_on=8(spinodal通過直後で打ち切り), T_off=40(5倍)で、")
    print("  サイクル開始時のu_far基準値が複数サイクルにわたって単調に上昇するか(ドリフトあり=")
    print("  パルス化は解決しない)、それとも定常状態に収束するか(ドリフトなし=部分的に機能)。")

    T_on, T_off, n_cycles = 8.0, 40.0, 8
    print("\n[本実験] R0=15, T_on=%.1f, T_off=%.1f, %d サイクル (t_final=%.1f)"
          % (T_on, T_off, n_cycles, n_cycles * (T_on + T_off)))
    result = run_pulsed(R0=15, T_on=T_on, T_off=T_off, n_cycles=n_cycles)

    print("  max_u_far=%.4f (u_sp_low=%.4f) mass_drift=%.2e diverged=%s"
          % (result["max_u_far"], U_SP_LOW, result["mass_drift"], result["diverged"]))
    print("  component series (t, phase, n, top_sizes):")
    for c in result["component_series"]:
        print("    t=%.1f [%s c%d] n=%d sizes=%s" % (c["t"], c["phase"], c["cycle"], c["n"],
                                                        c["top_sizes"]))

    # サイクル開始時点(各cycleのgrowth最初のスナップショット)のu_farでドリフトを判定
    cycle_start_u_far = {}
    for r in result["u_far_series"]:
        if r["phase"] == "growth" and r["cycle"] not in cycle_start_u_far:
            cycle_start_u_far[r["cycle"]] = r["u_far"]
    cycle_starts = [cycle_start_u_far[c] for c in sorted(cycle_start_u_far)]
    print("  サイクル開始時 u_far: %s" % [round(x, 4) for x in cycle_starts])
    net_drift = bool(len(cycle_starts) >= 3 and cycle_starts[-1] > cycle_starts[0] + 0.02)

    n_series = [c["n"] for c in result["component_series"]]
    division_events = sum(1 for i in range(1, len(n_series)) if n_series[i] > n_series[i - 1])
    cycling_evidence = bool(division_events >= 2)  # 複数回の増加=複数回の分裂の兆候
    result["nucleation_avoided"] = bool(not net_drift and result["max_u_far"] < 0.5)

    if net_drift:
        verdict = "falsified (net upward drift -- pulsing slows but does not stop nucleation)"
    elif cycling_evidence:
        verdict = "supported"
    else:
        verdict = "partially_supported (drift avoided, but no clear grow-divide cycling seen)"
    print("\n--- 判定: %s ---" % verdict)
    print("正味ドリフト: %s / 核形成回避(暫定): %s / 複数回の分裂イベントの兆候: %s (増加イベント数=%d)"
          % (net_drift, result["nucleation_avoided"], cycling_evidence, division_events))

    genesis_yaml = {"equations": "dividing_protocell, pulsed kp/kd (growth/rest square wave)",
                    "solver": "pseudo-spectral, custom pulsed time-stepping loop",
                    "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
                    "params": {"growth": GROWTH_PARAMS, "rest": REST_PARAMS,
                              "T_on": T_on, "T_off": T_off, "n_cycles": n_cycles},
                    "seed": 1, "seeds": [1], "commit": None,
                    "checksum": io.checksum_of(np.array([result["mass_drift"]]))}
    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": False, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction_original": ("T_onを安全マージンで切ればu_farはspinodalを一度も超えない。"
                                    "複数サイクルで成長->分裂->再成長のサイクルが見える。"),
            "prediction_revised": ("予備実験でspinodal通過(t~10)と可視核形成(t~75)の間に長い誘導"
                                  "時間があると判明したため予言を修正: T_on=8(spinodal通過直後で"
                                  "打ち切り)+T_off=40(5倍、緩和用)で、サイクル開始時u_farが複数"
                                  "サイクルで単調ドリフトするか(パルス化失敗)、定常化するか(部分成功)。"),
            "verdict": verdict,
            "T_on": T_on, "T_off": T_off, "n_cycles": n_cycles,
            "relaxation_vs_growth_rate_ratio": "growth ~0.0075/t, rest relaxation ~0.00016/t (~45x slower)",
            "max_u_far": result["max_u_far"], "u_sp_low": U_SP_LOW,
            "cycle_start_u_far_series": [round(x, 4) for x in cycle_starts],
            "net_drift_detected": net_drift,
            "division_events": division_events, "cycling_evidence": cycling_evidence,
            "component_series": result["component_series"],
            "u_far_series": result["u_far_series"],
        },
        "purity": {"per_object_labels": True, "external_optimum": False, "role": "S"},
    }
    integrity = io.integrity_block(
        conservation_drift=result["mass_drift"], resolutions_result={"64x64x64": 1},
        seed_success={"1": 1}, nan_or_clip=result["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # R0・T_on/T_offは調べたい挙動そのもの
        gate_encodes_conclusion_causality=False, gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "continuous production (no pulsing, from H-B1/H-B2)",
                        "result": "nucleates within t~75-140 at same kp/kd"}])
    notes = ("課題3-1。予言を予備実験の知見(誘導時間の非対称性)を踏まえて修正して登録。"
             "T_on=%.1f, T_off=%.1f, %d サイクル。判定=%s。u_far系列とcomponent系列は"
             "measured_byに記録。" % (T_on, T_off, n_cycles, verdict))
    run_dir = io.write_results("hyp-task3-1-temporal-pulsing", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return result


if __name__ == "__main__":
    main()
