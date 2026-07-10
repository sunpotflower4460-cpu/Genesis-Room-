# NOTES — g001-cgl-3d-b0c0-long-seed0001

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
    "variance_growth": 0.010619,
    "structure_factor_peak_k": 1.0,
    "structure_factor_prominence": 41.99999,
    "correlation_length": 1.0,
    "defect_count": 0,
    "defect_count_series": [
      48697,
      844,
      455,
      294,
      168,
      71,
      24,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0
    ],
    "defect_tail_mean": 0.0,
    "defect_dynamics_still_coarsening": false,
    "defect_dynamics_extinct": true,
    "defect_dynamics_trend_relative": 0.0,
    "defect_dynamics_residual_relative_std": 0.0
  }
}
```

## 出なかったもの・未達
reached_level=1, candidate_level=2

## EXPLORATION 所見
依頼2 [検証/探索]: (b,c)={'b': 0.0, 'c': 0.0}, t_final=400 (steps=8000, dt=0.0500), grid=(48, 48, 48)。3D の渦線再結合・結び目/絡み合いは common/diagnostics.winding_defect_count の貫通面カウント（線本数の下限プロキシ、docstring 記載）のみで評価しており、明示的な線追跡・結び目検出は未実装（未探索、依頼書の「足したら明記」に対応）。[検証] 標準 t_final=100 の run では tail がまだ単調減少中で「凍結」と断定できなかったため、t_final=400 まで延長して本当に定常/凍結に達するか確認する。
