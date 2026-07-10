# NOTES — g001-cgl-3d-b1p0-c0p5-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": false,
    "spontaneous_motion": false,
    "circulation": false,
    "persistent_individuality": false,
    "co_differentiation": false,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "variance_growth": 0.013074,
    "structure_factor_peak_k": 1.0,
    "structure_factor_prominence": 41.848029,
    "correlation_length": 1.0,
    "defect_count": 0,
    "defect_count_series": [
      53100,
      0,
      824,
      849,
      634,
      543,
      466,
      395,
      337,
      279,
      227,
      182,
      132,
      89,
      76,
      56,
      32,
      0,
      0,
      0,
      0,
      0
    ],
    "defect_tail_mean": 51.545,
    "defect_dynamics_still_coarsening": true,
    "defect_dynamics_extinct": false,
    "defect_dynamics_trend_relative": 3.7113,
    "defect_dynamics_residual_relative_std": 0.4497
  }
}
```

## 出なかったもの・未達
reached_level=1, candidate_level=2

## EXPLORATION 所見
依頼2 [検証/探索]: (b,c)={'b': 1.0, 'c': 0.5}, t_final=100 (steps=3011, dt=0.0332), grid=(48, 48, 48)。3D の渦線再結合・結び目/絡み合いは common/diagnostics.winding_defect_count の貫通面カウント（線本数の下限プロキシ、docstring 記載）のみで評価しており、明示的な線追跡・結び目検出は未実装（未探索、依頼書の「足したら明記」に対応）。[探索] (b,c)=(1.0,0.5), 1+bc=1.50 (安定域)。
