#!/usr/bin/env python3
"""genesis/ic_search_followup.py -- 大量ランダムIC探索の二段構え(2)＋law_variant余力枠(3)。

(2) 上位候補の深い評価: 台帳の上位(score_max順)を、より大きい格子・より長いt_finalで再評価し、
    「多個体構造が有限サイズによる一時的アーチファクトか、それとも本当に持続しないのか」を
    確かめる(finite-size hypothesis test)。多個体(blobs in [2,12])が持続した時間の割合を
    測り、compute_level_reportでreached_levelも算出する。

(3) law_variant sweep（余力枠、より厳しい監査つき）: 法則(b,c)を変えるのはIC変更よりも重い
    決定であり、docs/AI_EXPERIMENT_POLICY.md §4の追加監査を明示する:
      - 物理的由来: CGLはHopf分岐近傍の振動不安定性に対する標準正規形（既存g001_cgl_3d.pyと
        同じ方程式・同じ無次元化、係数の値だけを変える）。
      - 対称性: 任意の(b,c)でU(1)位相対称性・並進対称性は保たれる（変更しない）。
      - 次元整合性: b,cは無次元係数（既存ファイルと同じ無次元形、次元は変わらない）。
      - 保存則: CGLは散逸系で保存則を持たない（既存g001_cgl_3d.pyのintegrity_block同様、
        conservation_drift=0.0として記録、"該当なし"を明記）。
      - 既知極限: b=c=0で最も単純な緩和的極限(Allen-Cahn的)に帰着する既存の検証済みケース。
      - 非符号化(第8監査): (b,c)は宣言済みレンジからランダムに引く（結論を知って選ばない）。
        IC(family+params)は主探索の上位候補に固定し、法則だけを変える設計によって、
        「法則Xがスコアを上げる」という主張の決定的対照＝元の法則(b,c)=(1,-0.7)基準点、を
        常に揃えて比較する。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import emergence, io  # noqa: E402
from common.emergence_score import emergence_score  # noqa: E402
from genesis import g001_cgl_3d as cgl  # noqa: E402
from genesis import ic_families_3d as icf  # noqa: E402
from genesis.search_random_ic_3d import LAW_PARAMS_DEFAULT, run_cgl_with_ic  # noqa: E402


def deep_eval_candidate(family, params, big_shape, t_final, law_params=None, rng_seed=12345,
                         n_snapshots=25):
    """台帳のfamily+paramsを、より大きい格子・より長いt_finalで再評価する。ノイズの実現値
    (rng)は元のtrialと同一である必要はない(格子サイズが違えば同一ノイズ配列は原理的に
    不可能)——同じfamily+paramsという"IC生成レシピ"の頑健性を、独立ノイズ・大きい格子で
    確認するのが目的（同一seedでの単純な再現ではなく、頑健性チェック）。"""
    law_params = law_params or LAW_PARAMS_DEFAULT
    rng = np.random.default_rng(rng_seed)
    A0 = icf.GENERATORS[family](big_shape, rng, **params)
    snapshots, phys = run_cgl_with_ic(A0, t_final, law_params, n_snapshots=n_snapshots)

    per_step = []
    for i in range(1, len(snapshots)):
        s, d = emergence_score([snapshots[i - 1]["field"], snapshots[i]["field"]])
        per_step.append({"t": round(snapshots[i]["step"] * phys["dt"], 3), "score": s, "detail": d})
    best = max(per_step, key=lambda r: r["score"]) if per_step else {"t": 0.0, "score": 0.0, "detail": {}}
    blobs_series = [(r["t"], r["detail"].get("blobs", 0)) for r in per_step]
    sustained_multi = [1 for (_, b) in blobs_series if 2 <= b <= 12]
    sustained_frac = len(sustained_multi) / max(1, len(blobs_series))
    level_report = emergence.compute_level_report(snapshots, kind="ic_search_deep_eval")

    return {
        "family": family, "params": params, "grid": list(big_shape), "t_final": t_final,
        "law_params": law_params, "phys": phys, "score_max": best["score"], "score_max_t": best["t"],
        "score_max_detail": best["detail"], "blobs_series": blobs_series,
        "sustained_multi_grain_fraction": round(sustained_frac, 3),
        "reached_level": level_report["reached_level"], "detected": level_report["detected"],
        "final_field": snapshots[-1]["field"],
    }


def sample_bc(rng, b_range=(0.0, 1.5), c_range=(-2.0, 1.0)):
    return {"b": float(rng.uniform(*b_range)), "c": float(rng.uniform(*c_range))}


def law_variant_sweep(family, params, shape, t_final, n_variants=30, master_seed=777,
                       n_snapshots=13, baseline=None):
    """ICを固定(主探索の上位候補)し、(b,c)だけをランダムに振る。決定的対照として元のbaseline
    (b,c)=(1,-0.7)を必ず含め、常に同じICで比較する。"""
    baseline = baseline or LAW_PARAMS_DEFAULT
    rng = np.random.default_rng(master_seed)
    variants = [dict(baseline)] + [sample_bc(rng) for _ in range(n_variants)]
    results = []
    for bc in variants:
        ic_rng = np.random.default_rng(master_seed)  # 同じICノイズ実現(法則だけを変える対照設計)
        A0 = icf.GENERATORS[family](shape, ic_rng, **params)
        snapshots, phys = run_cgl_with_ic(A0, t_final, bc, n_snapshots=n_snapshots)
        per_step = []
        for i in range(1, len(snapshots)):
            s, d = emergence_score([snapshots[i - 1]["field"], snapshots[i]["field"]])
            per_step.append({"t": snapshots[i]["step"] * phys["dt"], "score": s, "detail": d})
        best = max(per_step, key=lambda r: r["score"]) if per_step else {"t": 0.0, "score": 0.0, "detail": {}}
        results.append({"b": bc["b"], "c": bc["c"], "bf": round(1 + bc["b"] * bc["c"], 4),
                        "score_max": best["score"], "score_max_t": round(best["t"], 2),
                        "diverged": phys["diverged"],
                        "is_baseline": bool(bc["b"] == baseline["b"] and bc["c"] == baseline["c"])})
    return results


if __name__ == "__main__":
    print("--- deep_eval_candidate smoke test ---")
    out = deep_eval_candidate("white_amp", {"amp": 0.44}, (40, 40, 40), 20.0)
    print("score_max=%.2f at t=%.1f, sustained_multi_grain_fraction=%.2f, reached_level=%d"
          % (out["score_max"], out["score_max_t"], out["sustained_multi_grain_fraction"],
             out["reached_level"]))

    print("\n--- law_variant_sweep smoke test (n=5) ---")
    res = law_variant_sweep("white_amp", {"amp": 0.44}, (24, 24, 24), 12.0, n_variants=5)
    for r in sorted(res, key=lambda r: r["score_max"], reverse=True):
        print("  b=%.2f c=%.2f (1+bc=%.2f) score_max=%.2f%s"
              % (r["b"], r["c"], r["bf"], r["score_max"], " [baseline]" if r["is_baseline"] else ""))
