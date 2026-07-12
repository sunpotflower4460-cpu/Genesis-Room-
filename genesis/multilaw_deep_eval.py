#!/usr/bin/env python3
"""genesis/multilaw_deep_eval.py -- 法則クラスまたぎ探索の二段構え深い評価。

CGLの教訓（短い観測窓での"持続"はFALSIFIEDだった、7+8監査で確認・訂正済み）を踏まえ、
Gray-Scott/Model Hの高速スクリーンでの"持続"候補も、同じ懐疑をもって長時間・大きい格子で
再評価する。CGLのdeep_eval_candidate(genesis/ic_search_followup.py)と同じ設計:
sustained_multi_grain_fraction(ここではsustained_multi_spot_fraction)、finite-sizeチェック。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import emergence  # noqa: E402
from common.emergence_score import emergence_score  # noqa: E402
from genesis.law_adapters import ADAPTERS  # noqa: E402


def deep_eval_candidate(law_name, family, params, big_shape, steps_total, rng_seed=999,
                        n_snapshots=25):
    adapter = ADAPTERS[law_name]
    rng = np.random.default_rng(rng_seed)
    state, used_params = adapter.ignite(big_shape, family, params, rng)
    k2_cache = adapter.k2_of(big_shape)
    snap_every = max(1, steps_total // (n_snapshots - 1))

    fields = [adapter.score_field(state).copy()]
    steps_list = [0]
    diverged = False
    done = 0
    while done < steps_total:
        chunk = min(snap_every, steps_total - done)
        state, diverged = adapter.step_n(state, adapter.default_dt, adapter.default_law_params,
                                          chunk, k2_cache=k2_cache)
        done += chunk
        fields.append(adapter.score_field(state).copy())
        steps_list.append(done)
        if diverged:
            break

    per_step = []
    for i in range(1, len(fields)):
        s, d = emergence_score([fields[i - 1], fields[i]])
        per_step.append({"t": steps_list[i], "score": s, "detail": d})
    best = max(per_step, key=lambda r: r["score"]) if per_step else {"t": 0, "score": 0.0, "detail": {}}
    blobs_series = [(r["t"], r["detail"].get("blobs", 0)) for r in per_step]
    sustained = [1 for (_, b) in blobs_series if 2 <= b <= 12]
    sustained_frac = len(sustained) / max(1, len(blobs_series))

    field_snapshots = [{"step": steps_list[i], "field": fields[i]} for i in range(len(fields))]
    level_report = emergence.compute_level_report(field_snapshots, kind="multilaw_deep_eval")

    return {
        "law": law_name, "family": family, "params": used_params, "grid": list(big_shape),
        "steps_total": done, "diverged": diverged, "score_max": best["score"],
        "score_max_t": best["t"], "blobs_series": blobs_series,
        "sustained_multi_spot_fraction": round(sustained_frac, 3),
        "reached_level": level_report["reached_level"], "detected": level_report["detected"],
        "final_field": fields[-1],
    }


if __name__ == "__main__":
    out = deep_eval_candidate("gray_scott", "patches", {"n_patches": 2, "radius_lo": 2.93,
                              "radius_hi": 3.29, "amp": 0.68}, (48, 48, 48), 8000, rng_seed=999)
    print("score_max=%.2f(t=%d) sustained_frac=%.2f reached_level=%d"
          % (out["score_max"], out["score_max_t"], out["sustained_multi_spot_fraction"],
             out["reached_level"]))
    print("blobs_series:", out["blobs_series"])
