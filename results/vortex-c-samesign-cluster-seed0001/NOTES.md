# NOTES — vortex-c-samesign-cluster-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": true,
    "spontaneous_motion": false,
    "circulation": true,
    "persistent_individuality": true,
    "co_differentiation": true,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "config": "samesign_cluster_2D",
    "n_vortex_points": 7,
    "bond_len": 6.93,
    "n_2cells": 6,
    "triangles_per_vertex": 0.857,
    "correlation_dim": 2.99,
    "ball_dim": null,
    "euler_char": 1,
    "n_components": 1,
    "n_edges": 12,
    "dimension_lifted_1_to_2": true,
    "ball_growth": [
      [
        0,
        1.0
      ],
      [
        1,
        4.429
      ],
      [
        2,
        7.0
      ],
      [
        3,
        7.0
      ],
      [
        4,
        7.0
      ],
      [
        5,
        7.0
      ],
      [
        6,
        7.0
      ]
    ],
    "energy_drift_frac": -1e-05,
    "honest_floor": "correlation_dim=位置(埋め込み)の幾何次元(固定格子); n_2cells=関係が三角形を閉じる創発量。空間自体の創発ではない。"
  }
}
```

## 出なかったもの・未達
reached_level=5, candidate_level=6

## EXPLORATION 所見
渦-依頼C [samesign_cluster_2D] t_final=20, grid=(64, 64, 64)。n_points=7, bond_len=6.93。n_2cells=6, corr_dim=2.99, ball_dim=None, components=1, dimension_lifted=True, E_drift=-0.0000。正例2: 共生クラスタ(7渦)の関係が6三角形を閉じ次元->2の予言。
