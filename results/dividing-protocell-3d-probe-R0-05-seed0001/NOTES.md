# NOTES — dividing-protocell-3d-probe-R0-05-seed0001

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
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "mass_drift": -0.0,
    "R0": 5,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 23,
    "division_time": 115.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 66,
    "significant_max_count": 66,
    "raw_final_count": 89,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      1664.0,
      1662.0,
      1660.0,
      1257.0,
      1255.0,
      1250.0,
      1250.0,
      1247.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=5.0 probe, t_final=150, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 1), (85, 1), (90, 1), (95, 1), (100, 1), (105, 1), (110, 2), (115, 19), (120, 72), (125, 86), (130, 95), (135, 87), (140, 86), (145, 87), (150, 89)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 1), (85, 1), (90, 1), (95, 1), (100, 1), (105, 1), (110, 1), (115, 7), (120, 8), (125, 51), (130, 34), (135, 34), (140, 57), (145, 63), (150, 66)]。division_mode=stable_no_division, daughter_size_ratio=None。R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs 多断片化）を特定するための決定的対照。
