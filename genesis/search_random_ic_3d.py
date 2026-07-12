#!/usr/bin/env python3
"""genesis/search_random_ic_3d.py -- 「これは物理だ」を保ったまま大量ランダムIC探索(3D)。

非交渉の前提（変えない）:
  - 法則は固定: genesis/g001_cgl_3d.py の complex Ginzburg-Landau, step()/stable_dt() をそのまま
    使う。(b,c) は既定でCLAUDEの2D予備調査と同じ検証済みregime (b=1.0, c=-0.7) に固定する
    (IC-onlyフェーズ)。law_variant は別関数(law_variant_sweep)でのみ、より厳しい監査つきで行う。
  - 変えるのは初期条件だけ: genesis/ic_families_3d.py の7 family からランダムサンプル。
  - t=0から放置: run中に手を入れない、完成形を置かない(各familyのdocstring参照、第8監査)。
  - スコアは common/emergence_score.py (2D検証済みを3D一般化) ＋ common/emergence.py の
    reached_level(Level 0-3を機械判定)を併用、台帳に両方記録。
  - 全trialを台帳(search_ledger/*.jsonl, append-only)に記録。失敗(死んだ/爆発/低スコア)も。

再現性: 単一のマスターRNG(np.random.default_rng(master_seed))から全ての乱択(family選択・
パラメータ・IC生成ノイズ)を順番に引く。同じmaster_seed・同じnなら、bandit重み付けを含めて
全trialが決定的に再現される(重みはtrialごとのオンライン更新のみに依存し、非決定的な外部
状態を持たない)。
"""

import argparse
import json
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import diagnostics as diag, emergence, io  # noqa: E402
from common.emergence_score import emergence_score  # noqa: E402
from genesis import g001_cgl_3d as cgl  # noqa: E402
from genesis import ic_families_3d as icf  # noqa: E402

LEDGER_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "search_ledger"))
LAW_PARAMS_DEFAULT = {"b": 1.0, "c": -0.7}  # Claude 2D予備調査の検証済みregime、法則固定


class Bandit:
    """epsilon-greedy softmax bandit（family選択）。事前バイアスで予備調査の上位
    (spectral_powerlaw, bandpass)を初期優遇しつつ、3Dでの実測平均スコアに応じて重みを
    オンライン更新する(3Dでランクが変わりうる、というhonest floorに対応)。"""

    def __init__(self, families, warmup_rounds=3, epsilon=0.25, temperature=1.5, prior_bias=None):
        self.families = list(families)
        self.counts = {f: 0 for f in self.families}
        self.sum_scores = {f: 0.0 for f in self.families}
        self.warmup_total = warmup_rounds * len(self.families)
        self.epsilon = epsilon
        self.temperature = temperature
        self.prior_bias = prior_bias or {}

    def mean(self, f):
        return self.sum_scores[f] / self.counts[f] if self.counts[f] > 0 else 0.0

    def choose(self, trial_idx, rng):
        if trial_idx < self.warmup_total:
            return self.families[trial_idx % len(self.families)]
        unsampled = [f for f in self.families if self.counts[f] == 0]
        if unsampled:
            return str(rng.choice(unsampled))
        if rng.random() < self.epsilon:
            return str(rng.choice(self.families))
        means = np.array([self.mean(f) for f in self.families])
        biases = np.array([self.prior_bias.get(f, 1.0) for f in self.families])
        scores = means * biases
        scores = scores - scores.max()
        exp = np.exp(scores / self.temperature)
        probs = exp / exp.sum()
        return str(rng.choice(self.families, p=probs))

    def update(self, family, score):
        self.counts[family] += 1
        self.sum_scores[family] += score

    def snapshot(self):
        return {f: {"n": self.counts[f], "mean_score": round(self.mean(f), 3)} for f in self.families}


def run_cgl_with_ic(A0, t_final, law_params, dt=None, n_snapshots=4):
    """法則固定(G001 CGL)でA0からt=0放置。cgl.step/stable_dtをそのまま使う(法則を変えない)。
    n_snapshots枚(等間隔、t=0と最終を含む)だけ保持してメモリを節約する。"""
    ndim = A0.ndim
    dt = dt if dt is not None else cgl.stable_dt(ndim, law_params["b"])
    steps = int(round(t_final / dt))
    snap_every = max(1, steps // (n_snapshots - 1)) if n_snapshots > 1 else steps
    A = A0.copy()
    snapshots = [{"step": 0, "field": A.copy()}]
    diverged = False
    for t in range(steps):
        A = cgl.step(A, dt, law_params)
        if not np.all(np.isfinite(A)):
            diverged = True
            break
        if (t + 1) % snap_every == 0 or t == steps - 1:
            snapshots.append({"step": t + 1, "field": A.copy()})
    phys = {"diverged": diverged, "dt": dt, "steps": steps}
    return snapshots, phys


def run_one_trial(family, shape, law_params, t_final, rng, n_snapshots=13):
    """スコアは軌道全体の"最良の瞬間"(score_max)でランクする——最終状態(score_final)だけを
    見ると、短時間の多粒状構造が速い粗大化で単一の巨大粒(box全体)へ収束した後を見てしまい、
    「途中でどこまで登ったか」を見落とす(t_final選びに敏感な測定アーチファクト、予備調査で
    発見・訂正)。churn(変化率)は常に隣接する2枚のスナップショットだけで測る(間隔が広いと
    churnを過小評価し、まだ落ち着いていない状態を誤って"persistent"と判定してしまうバグを
    避けるため)。reached_level_screen は最良の瞬間までの前方部分軌道から計算する。"""
    A0, params = icf.make_ic(family, shape, rng)
    snapshots, phys = run_cgl_with_ic(A0, t_final, law_params, n_snapshots=n_snapshots)
    per_step = []
    for i in range(1, len(snapshots)):
        s, d = emergence_score([snapshots[i - 1]["field"], snapshots[i]["field"]])
        per_step.append({"t": snapshots[i]["step"] * phys["dt"], "score": s, "detail": d, "idx": i})
    if not per_step:
        per_step = [{"t": 0.0, "score": 0.0, "detail": {}, "idx": 0}]
    best = max(per_step, key=lambda r: r["score"])
    final = per_step[-1]
    level_report = emergence.compute_level_report(snapshots[:best["idx"] + 1], kind="ic_search_screen")
    return {
        "family": family, "params": params,
        "score": best["score"], "score_detail": best["detail"], "score_t": round(best["t"], 3),
        "score_final": final["score"], "score_final_detail": final["detail"],
        "reached_level_screen": level_report["reached_level"],
        "detected_screen": level_report["detected"],
        "per_step_series": [(round(r["t"], 2), r["score"], r["detail"].get("blobs", 0)) for r in per_step],
        "phys": phys, "final_field": snapshots[-1]["field"], "A0": A0,
    }


def append_ledger(path, record):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "a") as f:
        f.write(json.dumps(record) + "\n")


def _json_safe(obj):
    if isinstance(obj, dict):
        return {k: _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    if isinstance(obj, (np.floating,)):
        return float(obj)
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    return obj


def main_search(n=500, shape=(28, 28, 28), t_final=15.0, master_seed=0,
                 law_params=None, ledger_name=None, prior_bias=None, n_snapshots=13):
    law_params = law_params or LAW_PARAMS_DEFAULT
    ledger_name = ledger_name or ("random_ic_3d_trials_seed%04d.jsonl" % master_seed)
    ledger_path = os.path.join(LEDGER_DIR, ledger_name)
    rng = np.random.default_rng(master_seed)
    prior_bias = prior_bias or {"spectral_powerlaw": 1.5, "bandpass": 1.3}
    bandit = Bandit(icf.FAMILIES, prior_bias=prior_bias)

    results = []
    t0 = time.time()
    for i in range(n):
        family = bandit.choose(i, rng)
        out = run_one_trial(family, shape, law_params, t_final, rng, n_snapshots=n_snapshots)
        bandit.update(family, out["score"])
        checksum = io.checksum_of(out["final_field"])
        record = {
            "trial": i, "family": family, "params": out["params"],
            "grid": list(shape), "law_params": law_params, "t_final": t_final,
            "dt": out["phys"]["dt"], "steps": out["phys"]["steps"],
            "diverged": out["phys"]["diverged"],
            "score": out["score"], "score_detail": out["score_detail"], "score_t": out["score_t"],
            "score_final": out["score_final"], "score_final_detail": out["score_final_detail"],
            "reached_level_screen": out["reached_level_screen"],
            "detected_screen": out["detected_screen"],
            "per_step_series": out["per_step_series"], "checksum": checksum,
            "master_seed": master_seed,
        }
        append_ledger(ledger_path, _json_safe(record))
        results.append(record)
        if (i + 1) % max(1, n // 20) == 0:
            elapsed = time.time() - t0
            print("  [%d/%d] elapsed=%.1fs  bandit=%s" % (i + 1, n, elapsed, bandit.snapshot()))

    print("=== search done: %d trials in %.1fs ===" % (n, time.time() - t0))
    print("bandit final:", bandit.snapshot())
    return results, bandit, ledger_path


def rank_results(results, top_k=15):
    return sorted(results, key=lambda r: r["score"], reverse=True)[:top_k]


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=500)
    ap.add_argument("--grid", type=int, default=28)
    ap.add_argument("--t_final", type=float, default=15.0)
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    shape = (args.grid, args.grid, args.grid)
    print("=== 大量ランダムIC探索(3D) === grid=%s n=%d t_final=%.0f seed=%d law=%s"
          % (shape, args.n, args.t_final, args.seed, LAW_PARAMS_DEFAULT))
    results, bandit, ledger_path = main_search(n=args.n, shape=shape, t_final=args.t_final,
                                                master_seed=args.seed)
    print("\n=== top %d (ranked by score_max, the best moment along each trajectory) ===" % 15)
    for r in rank_results(results, 15):
        print("  trial=%d family=%-18s score=%.2f(t=%.1f) score_final=%.2f reached=%d params=%s"
              % (r["trial"], r["family"], r["score"], r["score_t"], r["score_final"],
                 r["reached_level_screen"], r["params"]))
    print("\nledger:", ledger_path)
