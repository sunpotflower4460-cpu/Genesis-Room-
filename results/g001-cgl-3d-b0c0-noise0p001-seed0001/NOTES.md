# NOTES — g001-cgl-3d-b0c0-noise0p001-seed0001

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
    "variance_growth": 0.010426,
    "structure_factor_peak_k": 1.0,
    "structure_factor_prominence": 39.448193,
    "correlation_length": 1.0,
    "defect_count": 109,
    "defect_count_series": [
      48678,
      0,
      0,
      969,
      778,
      642,
      523,
      476,
      427,
      385,
      364,
      345,
      312,
      292,
      277,
      232,
      198,
      148,
      136,
      122,
      109
    ],
    "defect_tail_mean": 230.455,
    "defect_dynamics_still_coarsening": true,
    "defect_dynamics_extinct": false,
    "defect_dynamics_trend_relative": 1.3287,
    "defect_dynamics_residual_relative_std": 0.0526
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
依頼2 [検証/探索]: (b,c)={'b': 0.0, 'c': 0.0}, t_final=100 (steps=2000, dt=0.0500), grid=(48, 48, 48)。3D の渦線再結合・結び目/絡み合いは common/diagnostics.winding_defect_count の貫通面カウント（線本数の下限プロキシ、docstring 記載）のみで評価しており、明示的な線追跡・結び目検出は未実装（未探索、依頼書の「足したら明記」に対応）。[探索] noise_amplitude=0.001（既定 0.01 との比較、KZ 的な初期条件依存の代理チェック。真のクエンチ速度スキャンではない点に注意）。
