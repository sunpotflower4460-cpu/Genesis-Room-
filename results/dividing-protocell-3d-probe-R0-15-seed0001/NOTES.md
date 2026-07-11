# NOTES — dividing-protocell-3d-probe-R0-15-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": false,
    "localization": true,
    "spontaneous_motion": false,
    "circulation": false,
    "persistent_individuality": true,
    "co_differentiation": false,
    "self_maintaining_closure": false,
    "growth_division_inheritance": true,
    "selection_open_ended": false
  },
  "measured_by": {
    "mass_drift": 1e-15,
    "R0": 15,
    "division_mode": "multi_fragmentation",
    "division_snapshot_index": 15,
    "division_time": 75.0,
    "daughter_size_ratio": 0.4397,
    "significant_initial_count": 1,
    "significant_final_count": 20,
    "significant_max_count": 20,
    "raw_final_count": 109,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      2607.0,
      1314.0,
      1313.0,
      1312.0,
      801.0,
      800.0,
      796.0,
      791.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=7, candidate_level=8

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=15.0 probe, t_final=150, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 2), (20, 2), (25, 2), (30, 2), (35, 2), (40, 2), (45, 2), (50, 2), (55, 2), (60, 2), (65, 2), (70, 2), (75, 9), (80, 18), (85, 54), (90, 12), (95, 26), (100, 21), (105, 39), (110, 40), (115, 63), (120, 63), (125, 63), (130, 63), (135, 74), (140, 106), (145, 109), (150, 109)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 2), (80, 2), (85, 2), (90, 2), (95, 5), (100, 6), (105, 6), (110, 6), (115, 9), (120, 9), (125, 15), (130, 15), (135, 14), (140, 19), (145, 13), (150, 20)]。division_mode=multi_fragmentation, daughter_size_ratio=0.4397。R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs 多断片化）を特定するための決定的対照。
