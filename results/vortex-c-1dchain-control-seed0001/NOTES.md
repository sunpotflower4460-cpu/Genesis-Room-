# NOTES — vortex-c-1dchain-control-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": true,
    "spontaneous_motion": false,
    "circulation": true,
    "persistent_individuality": true,
    "co_differentiation": false,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "config": "chain_1D",
    "n_vortex_points": 4,
    "bond_len": 8.277,
    "n_2cells": 0,
    "triangles_per_vertex": 0.0,
    "correlation_dim": 1.042,
    "ball_dim": 0.485,
    "euler_char": 1,
    "n_components": 1,
    "n_edges": 3,
    "dimension_lifted_1_to_2": false,
    "ball_growth": [
      [
        0,
        1.0
      ],
      [
        1,
        2.5
      ],
      [
        2,
        3.5
      ],
      [
        3,
        4.0
      ],
      [
        4,
        4.0
      ],
      [
        5,
        4.0
      ],
      [
        6,
        4.0
      ]
    ],
    "energy_drift_frac": -6e-05,
    "honest_floor": "correlation_dim=位置(埋め込み)の幾何次元(固定格子); n_2cells=関係が三角形を閉じる創発量。空間自体の創発ではない。"
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
渦-依頼C [chain_1D] t_final=20, grid=(64, 64, 64)。n_points=4, bond_len=8.28。n_2cells=0, corr_dim=1.042, ball_dim=0.485, components=1, dimension_lifted=False, E_drift=-0.0001。決定的対照2: 関係が1次元の紐、2-cell立たずの予言。
