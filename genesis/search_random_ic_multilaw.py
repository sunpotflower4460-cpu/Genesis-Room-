#!/usr/bin/env python3
"""genesis/search_random_ic_multilaw.py -- 法則クラスをまたぐ大量ランダムIC探索。

CGL(genesis/search_random_ic_3d.py)の(b,c)law_variant調査で、Benjamin-Feir不安定は
"乱流で"粗大化を防ぐだけで清潔な持続個体を作らないことが分かった(7+8監査で指摘・訂正済み)。
うえきの選択で、法則を(b,c)調整でなく法則クラス自体(Gray-Scott反応拡散・Model H相分離+流れ)
を切り替えて同じ大量ランダムIC探索を回す。genesis/law_adapters.py の統一インタフェースを使い、
CGL版と同じ Bandit・スコア・台帳の設計をそのまま再利用する。
"""

import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import emergence, io  # noqa: E402
from common.emergence_score import emergence_score  # noqa: E402
from genesis import ic_families_real as icr  # noqa: E402
from genesis.law_adapters import ADAPTERS  # noqa: E402
from genesis.search_random_ic_3d import Bandit, append_ledger, _json_safe  # noqa: E402

LEDGER_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "search_ledger"))


def run_one_trial(adapter, family, shape, law_params, steps_total, rng, n_snapshots=13):
    state, params = adapter.ignite(shape, family, None, rng)
    k2_cache = adapter.k2_of(shape)
    snap_every = max(1, steps_total // (n_snapshots - 1))
    fields = [adapter.score_field(state).copy()]
    steps_list = [0]
    diverged = False
    done = 0
    while done < steps_total:
        chunk = min(snap_every, steps_total - done)
        state, diverged = adapter.step_n(state, adapter.default_dt, law_params, chunk, k2_cache=k2_cache)
        done += chunk
        fields.append(adapter.score_field(state).copy())
        steps_list.append(done)
        if diverged:
            break

    per_step = []
    for i in range(1, len(fields)):
        s, d = emergence_score([fields[i - 1], fields[i]])
        per_step.append({"t": steps_list[i], "score": s, "detail": d})
    if not per_step:
        per_step = [{"t": 0, "score": 0.0, "detail": {}}]
    best = max(per_step, key=lambda r: r["score"])
    final = per_step[-1]

    field_snapshots = [{"step": steps_list[i], "field": fields[i]} for i in range(len(fields))]
    best_idx = next(i for i, r in enumerate(per_step) if r is best) + 1
    level_report = emergence.compute_level_report(field_snapshots[:best_idx + 1], kind="ic_search_multilaw")

    return {
        "family": family, "params": params, "score": best["score"], "score_detail": best["detail"],
        "score_t": best["t"], "score_final": final["score"], "score_final_detail": final["detail"],
        "reached_level_screen": level_report["reached_level"], "detected_screen": level_report["detected"],
        "per_step_series": [(r["t"], r["score"], r["detail"].get("blobs", 0)) for r in per_step],
        "diverged": diverged, "final_field": fields[-1], "steps_total": done,
    }


def main_search(law_name, n=150, shape=None, steps_total=2000, master_seed=0, prior_bias=None,
                 n_snapshots=13):
    adapter = ADAPTERS[law_name]
    shape = shape or adapter.default_grid
    law_params = adapter.default_law_params
    ledger_path = os.path.join(LEDGER_DIR, "random_ic_%s_trials_seed%04d.jsonl" % (law_name, master_seed))
    rng = np.random.default_rng(master_seed)
    bandit = Bandit(icr.FAMILIES, prior_bias=prior_bias or {})

    results = []
    t0 = time.time()
    for i in range(n):
        family = bandit.choose(i, rng)
        out = run_one_trial(adapter, family, shape, law_params, steps_total, rng, n_snapshots=n_snapshots)
        bandit.update(family, out["score"])
        checksum = io.checksum_of(out["final_field"])
        record = {
            "trial": i, "law": law_name, "family": family, "params": out["params"],
            "grid": list(shape), "law_params": law_params, "steps_total": out["steps_total"],
            "diverged": out["diverged"],
            "score": out["score"], "score_detail": out["score_detail"], "score_t": out["score_t"],
            "score_final": out["score_final"], "score_final_detail": out["score_final_detail"],
            "reached_level_screen": out["reached_level_screen"], "detected_screen": out["detected_screen"],
            "per_step_series": out["per_step_series"], "checksum": checksum, "master_seed": master_seed,
        }
        append_ledger(ledger_path, _json_safe(record))
        results.append(record)
        if (i + 1) % max(1, n // 15) == 0:
            elapsed = time.time() - t0
            print("  [%d/%d] elapsed=%.1fs bandit=%s" % (i + 1, n, elapsed, bandit.snapshot()))

    print("=== %s search done: %d trials in %.1fs ===" % (law_name, n, time.time() - t0))
    print("bandit final:", bandit.snapshot())
    return results, bandit, ledger_path


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--law", choices=["gray_scott", "model_h"], required=True)
    ap.add_argument("--n", type=int, default=150)
    ap.add_argument("--steps_total", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    print("=== 法則クラスをまたぐ探索: %s ===" % args.law)
    results, bandit, ledger_path = main_search(args.law, n=args.n, steps_total=args.steps_total,
                                                master_seed=args.seed)
    top = sorted(results, key=lambda r: -r["score"])[:15]
    print("\n=== top 15 ===")
    for r in top:
        print("  trial=%d family=%-20s score=%.2f(t=%d) score_final=%.2f reached=%d"
              % (r["trial"], r["family"], r["score"], r["score_t"], r["score_final"],
                 r["reached_level_screen"]))
    print("\nledger:", ledger_path)
