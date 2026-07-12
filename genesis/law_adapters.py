#!/usr/bin/env python3
"""genesis/law_adapters.py -- 法則クラスをまたぐ大量ランダムIC探索のための統一インタフェース。

CGL(依頼A/B/C/Dの前回ラウンド)の(b,c)調整は天井5.7で頭打ちだった——BF不安定は"乱流で"
粗大化を防ぐのであって"清潔な持続個体"を作らない、という発見を受け、うえきの選択で「法則
クラスをまたぐ探索」を実行する: 既存の genesis/gray_scott_3d.py（反応拡散、膜なし）と
genesis/g003_model_h_3d.py（Cahn-Hilliard相分離+Navier-Stokes流れ）に、同じ大量ランダムIC
探索harnessを適用する。

各 adapter は以下を提供する（統一契約）:
  ignite(shape, family, params, rng) -> background_state    法則固有の背景に family 摂動を
                                                              合成した初期状態を作る
  step_n(state, dt, law_params, n)  -> state                n ステップ進める(法則固定)
  score_field(state) -> ndarray                              emergence_score に渡すND実数配列
                                                              (法則の"order parameter")
  default_law_params, default_grid, default_dt               法則ごとの既定値

第8監査: どの adapter も、family 摂動を「背景への摂動」として合成するだけで、完成した
パターン(GSのスポット・Model Hの液滴)を初期条件に置かない。GSは既知のbistabilityにより
空間相関のある摂動(ic_families_real の spectral_powerlaw/bandpass/patches/seeds等)でないと
点火しないことを予備確認済み(white_ampは相関なしのため点火しにくいと予想され、これ自体が
familyの弁別に使える)——u=1-v の対で摂動を合成するのは標準的なGS規約であり、パターンの
形状を指定するものではない。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from genesis import g003_model_h_3d as m_h  # noqa: E402
from genesis import gray_scott_3d as gs  # noqa: E402
from genesis import ic_families_real as icr  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402


# GSは可興奮媒質(bistable)であり、CGLと違って「どんな摂動でも育つ」わけではない——予備
# 確認で幅radius~6-12(28^3ボックスの21-43%!)の広いパッチは局所的な基質(u)を広範囲に
# 枯渇させてしまい、反応が持続する前に緩和して消える(点火しない)ことが分かった。既存
# genesis/gray_scott_3d.pyのmake_seed_initialはseed_radius=4という「局所的に集中した」
# 慣行を使っている——それに合わせ、GS専用に幅/半径のレンジを狭める(第8監査: 形は
# family任せのまま、GSの物理的な点火条件に合わせた"スケール"だけを調整。docs/
# AI_EXPERIMENT_POLICY.md §2で明示的に許可されている「空間スケール」の調整)。
GS_PARAM_RANGES = {
    "spectral_powerlaw": {"beta": (0.5, 4.0), "amp": (0.2, 0.6)},
    "bandpass": {"klo": (1.0, 6.0), "kw": (1.0, 4.0), "amp": (0.2, 0.6)},
    "white_amp": {"amp": (0.2, 0.6)},
    "seeds_signed": {"n_seeds": (1, 5), "width": (1.5, 4.0), "amp": (0.3, 0.7)},
    "line_or_ring_patches": {"n_defects": (1, 2), "amp": (0.3, 0.7)},
    "patches": {"n_patches": (1, 4), "radius_lo": (1.5, 3.0), "radius_hi": (3.0, 5.0),
                "amp": (0.3, 0.7)},
    "seeds": {"n_seeds": (1, 5), "width": (1.5, 4.0), "amp": (0.3, 0.7)},
}


class GrayScottAdapter:
    name = "gray_scott"
    default_law_params = dict(gs.DEFAULTS)  # F=0.0367,k=0.0649 (mitosis域、法則固定)
    default_grid = (28, 28, 28)
    default_dt = gs.DT

    def ignite(self, shape, family, params, rng):
        if params is not None:
            p = params
        else:
            # GS専用レンジでサンプル(局所的な"点火"に合うスケール、familyの形自体は同じ関数)
            ranges = GS_PARAM_RANGES[family]
            p = {}
            for key, (lo, hi) in ranges.items():
                if key in ("n_seeds", "n_defects", "n_patches"):
                    p[key] = int(rng.integers(lo, hi + 1))
                else:
                    p[key] = float(rng.uniform(lo, hi))
        pert = icr.GENERATORS[family](shape, rng, **p)
        v = np.clip(np.abs(pert), 0.0, 1.0)  # v>=0のGS規約、摂動の絶対値を使う
        u = np.clip(1.0 - v, 0.0, 1.0)  # 標準的な"局所的に反応物が減る"対のimprint規約
        return {"u": u, "v": v}, p

    def step_n(self, state, dt, law_params, n, k2_cache=None):
        u, v = state["u"], state["v"]
        _, k2 = k2_cache if k2_cache is not None else k_grid(u.shape)
        diverged = False
        for _ in range(n):
            u, v = gs.step(u, v, dt, law_params, k2)
            if not (np.all(np.isfinite(u)) and np.all(np.isfinite(v))):
                diverged = True
                break
        return {"u": u, "v": v}, diverged

    def score_field(self, state):
        return state["v"]

    def k2_of(self, shape):
        return k_grid(shape)


class ModelHAdapter:
    name = "model_h"
    default_law_params = dict(m_h.DEFAULTS)  # M,kappa,nu,C 固定(法則固定)
    default_grid = (28, 28, 28)
    default_dt = 0.02  # g003_model_h_3d.py 本体の t_final/steps 慣行に準拠した控えめなdt

    def ignite(self, shape, family, params, rng):
        pert, used_params = icr.make_ic(family, shape, rng, params=params)
        # 一様background(mean=0) + family摂動をphiとする(g003本体のmake_initialと同じ
        # "完成した液滴/界面は置かない"規約、摂動の大きさだけfamilyに委ねる)
        phi = pert
        ndim = len(shape)
        u = [np.zeros(shape) for _ in range(ndim)]  # 流れは静止から(完成した流れ場を置かない)
        return {"phi": phi, "u": u}, used_params

    def step_n(self, state, dt, law_params, n, k2_cache=None):
        from genesis.solvers import dealias_mask
        phi, u = state["phi"], state["u"]
        shape = phi.shape
        if k2_cache is not None:
            kk, k2 = k2_cache
        else:
            kk, k2 = k_grid(shape)
        dealias = dealias_mask(shape)
        diverged = False
        for _ in range(n):
            phi, u, _mu = m_h.step(phi, u, dt, law_params, kk, k2, dealias)
            if not np.all(np.isfinite(phi)):
                diverged = True
                break
        return {"phi": phi, "u": u}, diverged

    def score_field(self, state):
        return state["phi"]

    def k2_of(self, shape):
        from genesis.solvers import k_grid as kg
        return kg(shape)


ADAPTERS = {"gray_scott": GrayScottAdapter(), "model_h": ModelHAdapter()}


if __name__ == "__main__":
    # 点火の予備確認: 各lawで各familyがscore_fieldに構造を作るか
    from common.emergence_score import emergence_score
    for law_name, adapter in ADAPTERS.items():
        print("=== %s ===" % law_name)
        shape = (24, 24, 24)
        for fam in icr.FAMILIES:
            rng = np.random.default_rng(1)
            state, params = adapter.ignite(shape, fam, None, rng)
            f0 = adapter.score_field(state).copy()
            k2c = adapter.k2_of(shape)
            state, diverged = adapter.step_n(state, adapter.default_dt, adapter.default_law_params,
                                              200, k2_cache=k2c)
            f1 = adapter.score_field(state)
            score, detail = emergence_score([f0, f1])
            print("  %-20s diverged=%s score=%.2f blobs=%d relvar=%.3f field_mean=%.4f"
                  % (fam, diverged, score, detail.get("blobs", -1), detail.get("rel_var", -1), f1.mean()))
