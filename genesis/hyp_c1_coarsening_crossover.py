#!/usr/bin/env python3
"""H-C1: G003 の粗大化は流体力学クロスオーバー前（依頼1-3の48^3・t_final=100が小さすぎた）。

背景: 依頼1-3で測った粗大化指数 n~0.25 は Siggia（粘性域 n=1）未満だった。指数は初期/小ドメイン
（拡散的 n~1/3）から晩期/大ドメイン（流体力学長 L_h を超えると n~2/3 慣性 or 1 粘性）へ
クロスオーバーし、依頼1-3の設定ではクロスオーバー前だった、という仮説。

--- 予言（実行前に登録） ---
グリッドを 64^3（依頼1-3の48^3より大きい）へ、時間を t_final=300（依頼1-3の100の3倍）へ
延ばして局所指数 n(t) = d(log L)/d(log t) をスライディングウィンドウで測ると：
  1. Model H（C=2、既定）: L の成長とともに n(t) が時間とともに上昇する（拡散的域から
     流体力学域へのクロスオーバー）。
  2. C=0（純粋 CH、決定的対照）: n(t) は上昇せず、拡散的 n~1/3 近辺に留まる。

falsification: Model H でも n(t) が時間とともに上昇しない（横ばいまたは低下）なら、この格子・
時間スケールでは流体力学的クロスオーバーの兆候がない、という結論になり H-C1 は否定される
（計算資源の制約で 96^3/128^3 まで届いていない点は限界として明記する）。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import g003_model_h_3d as g003  # noqa: E402

SHAPE = (64, 64, 64)
DT = 0.05
STEPS = 6000  # t_final = 300 (依頼1-3の3倍)
SNAPSHOT_EVERY = 60  # 100 snapshots -- n(t) のスライディングフィットに十分な時間分解能


def local_exponent_series(times, lengths, window=9):
    """log(L) vs log(t) のスライディングウィンドウ線形回帰で局所指数 n(t) を計算する。"""
    pts = [(t, l) for t, l in zip(times, lengths) if t > 0 and l > 0]
    if len(pts) < window:
        return [], []
    logt = np.log([p[0] for p in pts])
    logl = np.log([p[1] for p in pts])
    t_centers, n_vals = [], []
    half = window // 2
    for i in range(half, len(pts) - half):
        sl = slice(i - half, i + half + 1)
        n = float(np.polyfit(logt[sl], logl[sl], 1)[0])
        t_centers.append(pts[i][0])
        n_vals.append(n)
    return t_centers, n_vals


def run_and_measure(params, seed=1):
    snapshots, phys = g003.run(SHAPE, STEPS, DT, seed, params=params, snapshot_every=SNAPSHOT_EVERY)
    times = [s["step"] * DT for s in snapshots]
    lengths = [diag.coarsening_length(s["field"]) for s in snapshots]
    t_centers, n_vals = local_exponent_series(times, lengths)
    return {
        "times": times, "lengths": lengths, "n_t_centers": t_centers, "n_t_values": n_vals,
        "n_early": float(np.mean(n_vals[:3])) if len(n_vals) >= 3 else None,
        "n_late": float(np.mean(n_vals[-3:])) if len(n_vals) >= 3 else None,
        "mass_drift": phys["mass_drift"], "diverged": phys["diverged"],
        "free_energy_drift": phys["free_energy_drift"],
    }


def main():
    print("=== H-C1: G003粗大化クロスオーバー (grid=64^3, t_final=300) ===")
    print("予言: Model H(C=2)はn(t)が時間とともに上昇、C=0対照はn~1/3に留まる")

    print("[Model H, C=2 既定]")
    mh = run_and_measure(dict(g003.DEFAULTS))
    print("  n(t): %s" % [round(n, 3) for n in mh["n_t_values"]])
    print("  n_early=%.3f n_late=%.3f mass_drift=%.2e" % (mh["n_early"], mh["n_late"], mh["mass_drift"]))

    print("[C=0 決定的対照（純粋CH）]")
    c0 = run_and_measure(dict(g003.DEFAULTS, C=0.0))
    print("  n(t): %s" % [round(n, 3) for n in c0["n_t_values"]])
    print("  n_early=%.3f n_late=%.3f mass_drift=%.2e" % (c0["n_early"], c0["n_late"], c0["mass_drift"]))

    mh_rises = bool(mh["n_late"] is not None and mh["n_early"] is not None
                    and mh["n_late"] > mh["n_early"] + 0.03)
    c0_flat = bool(c0["n_late"] is not None and c0["n_early"] is not None
                   and abs(c0["n_late"] - c0["n_early"]) < 0.1)
    verdict = "supported" if (mh_rises and not (c0["n_late"] and c0["n_late"] > c0["n_early"] + 0.03)) else "falsified"

    print("\n--- 判定: %s ---" % verdict)
    print("Model H で n(t) 上昇 (n_early->n_late, +0.03超): %s (%.3f -> %.3f)"
          % (mh_rises, mh["n_early"], mh["n_late"]))
    print("C=0対照は横ばい: %s (%.3f -> %.3f)" % (c0_flat, c0["n_early"], c0["n_late"]))

    report = {
        "reached_level": None, "candidate_level": None,
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": {}, "measured_by": {
            "prediction": "Model H(C=2)はn(t)が時間とともに上昇（クロスオーバー）、C=0はn~1/3に留まる",
            "verdict": verdict,
            "model_h": {k: v for k, v in mh.items() if k not in ("times", "lengths")},
            "control_C0": {k: v for k, v in c0.items() if k not in ("times", "lengths")},
            "model_h_n_series": mh["n_t_values"], "control_C0_n_series": c0["n_t_values"],
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": "E"},
    }
    integrity = io.integrity_block(
        conservation_drift=max(abs(mh["mass_drift"]), abs(c0["mass_drift"])),
        resolutions_result={"64x64x64_t300": mh["n_late"]}, seed_success={"1": 1},
        nan_or_clip=mh["diverged"] or c0["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "C=0 (pure CH, no hydrodynamics)",
                        "result": "n_early=%.3f n_late=%.3f" % (c0["n_early"], c0["n_late"])}])
    genesis_yaml = {"equations": "Model H vs C=0, long-time coarsening exponent crossover test",
                    "solver": "pseudo-spectral (g003_model_h_3d.run)",
                    "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
                    "params": "C=2(Model H) vs C=0(control)", "seed": 1, "seeds": [1], "commit": None,
                    "checksum": io.checksum_of(np.array(mh["n_t_values"] + c0["n_t_values"]))}
    notes = ("H-C1検証。優先順位4位。64^3・t_final=300（依頼1-3比: grid 1.19x/time 3x大）。"
             "計算資源制約で依頼書提案の96^3/128^3には未到達（限界として明記）。判定=%s。"
             % verdict)
    run_dir = io.write_results("hyp-c1-coarsening-crossover", genesis_yaml, report, integrity,
                                input_vs_output, figures={}, notes=notes)
    print("wrote %s" % run_dir)
    return mh, c0


if __name__ == "__main__":
    main()
