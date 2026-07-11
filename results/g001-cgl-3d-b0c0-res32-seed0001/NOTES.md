# NOTES — g001-cgl-3d-b0c0-res32-seed0001

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
    "variance_growth": 0.011727,
    "structure_factor_peak_k": 1.0,
    "structure_factor_prominence": 27.999861,
    "correlation_length": 1.0,
    "defect_count": 0,
    "defect_count_series": [
      16570,
      0,
      296,
      262,
      129,
      55,
      11,
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
依頼2 [検証/探索]: (b,c)={'b': 0.0, 'c': 0.0}, t_final=100 (steps=2000, dt=0.0500), grid=(32, 32, 32)。3D の渦線再結合・結び目/絡み合いは common/diagnostics.winding_defect_count の貫通面カウント（線本数の下限プロキシ、docstring 記載）のみで評価しており、明示的な線追跡・結び目検出は未実装（未探索、依頼書の「足したら明記」に対応）。解像度収束チェック（32^3 vs 48^3、同一 spacing=1.0 -- 注意: N のみ変えているため物理箱サイズも変わっており、厳密な意味での「同一領域を細かくする」解像度収束ではなく有限サイズ依存チェックに近い。
