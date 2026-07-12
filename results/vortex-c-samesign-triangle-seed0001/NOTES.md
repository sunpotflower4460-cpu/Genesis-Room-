# NOTES — vortex-c-samesign-triangle-seed0001

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
    "config": "samesign_triangle",
    "n_vortex_points": 3,
    "bond_len": 9.899,
    "n_2cells": 1,
    "triangles_per_vertex": 0.333,
    "correlation_dim": null,
    "ball_dim": null,
    "euler_char": 1,
    "n_components": 1,
    "n_edges": 3,
    "dimension_lifted_1_to_2": true,
    "ball_growth": [
      [
        0,
        1.0
      ],
      [
        1,
        3.0
      ],
      [
        2,
        3.0
      ],
      [
        3,
        3.0
      ],
      [
        4,
        3.0
      ],
      [
        5,
        3.0
      ],
      [
        6,
        3.0
      ]
    ],
    "energy_drift_frac": -2e-05,
    "honest_floor": "correlation_dim=位置(埋め込み)の幾何次元(固定格子); n_2cells=関係が三角形を閉じる創発量。空間自体の創発ではない。"
  }
}
```

## 出なかったもの・未達
reached_level=5, candidate_level=6

## EXPLORATION 所見
渦-依頼C [samesign_triangle] t_final=20, grid=(64, 64, 64)。n_points=3, bond_len=9.90。n_2cells=1, corr_dim=None, ball_dim=None, components=1, dimension_lifted=True, E_drift=-0.0000。正例1(最小): 共生する3渦の関係が三角形(1 2-cell)を閉じる予言。
