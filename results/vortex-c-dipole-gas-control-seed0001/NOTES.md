# NOTES — vortex-c-dipole-gas-control-seed0001

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
    "config": "dipole_gas",
    "n_vortex_points": 18,
    "bond_len": 4.427,
    "n_2cells": 0,
    "triangles_per_vertex": 0.0,
    "correlation_dim": 1.545,
    "ball_dim": 0.396,
    "euler_char": 7,
    "n_components": 7,
    "n_edges": 11,
    "dimension_lifted_1_to_2": false,
    "ball_growth": [
      [
        0,
        1.0
      ],
      [
        1,
        2.222
      ],
      [
        2,
        2.889
      ],
      [
        3,
        3.444
      ],
      [
        4,
        3.889
      ],
      [
        5,
        4.222
      ],
      [
        6,
        4.444
      ]
    ],
    "energy_drift_frac": -7e-05,
    "honest_floor": "correlation_dim=位置(埋め込み)の幾何次元(固定格子); n_2cells=関係が三角形を閉じる創発量。空間自体の創発ではない。"
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
渦-依頼C [dipole_gas] t_final=20, grid=(64, 64, 64)。n_points=18, bond_len=4.43。n_2cells=0, corr_dim=1.545, ball_dim=0.396, components=7, dimension_lifted=False, E_drift=-0.0001。決定的対照1: 関係が辺だけ、2-cellゼロの予言。
