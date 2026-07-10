# NOTES — g001-cgl-3d-b0p3-c0p3-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": true,
    "spontaneous_motion": false,
    "circulation": false,
    "persistent_individuality": false,
    "co_differentiation": false,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "variance_growth": 0.013345,
    "structure_factor_peak_k": 1.0,
    "structure_factor_prominence": 40.111644,
    "correlation_length": 1.0,
    "defect_count": 114,
    "defect_count_series": [
      48683,
      0,
      1076,
      998,
      810,
      662,
      565,
      502,
      476,
      446,
      429,
      408,
      384,
      367,
      355,
      331,
      312,
      265,
      215,
      163,
      114
    ],
    "defect_tail_mean": 303.909,
    "defect_dynamics_still_coarsening": true,
    "defect_dynamics_extinct": false,
    "defect_dynamics_trend_relative": 1.0888,
    "defect_dynamics_residual_relative_std": 0.0774
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
依頼2 [検証/探索]: (b,c)={'b': 0.3, 'c': 0.3}, t_final=100 (steps=2000, dt=0.0500), grid=(48, 48, 48)。3D の渦線再結合・結び目/絡み合いは common/diagnostics.winding_defect_count の貫通面カウント（線本数の下限プロキシ、docstring 記載）のみで評価しており、明示的な線追跡・結び目検出は未実装（未探索、依頼書の「足したら明記」に対応）。[探索] (b,c)=(0.3,0.3), 1+bc=1.09 (安定域)。
