# NOTES — dividing-protocell-3d-probe-R0-11-seed0001

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
    "mass_drift": -3e-15,
    "R0": 11,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 18,
    "division_time": 90.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 14,
    "significant_max_count": 19,
    "raw_final_count": 66,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      11748.0,
      848.0,
      843.0,
      838.0,
      837.0,
      835.0,
      829.0,
      826.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=11.0 probe, t_final=150, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 1), (85, 59), (90, 2), (95, 5), (100, 14), (105, 14), (110, 26), (115, 77), (120, 71), (125, 79), (130, 79), (135, 75), (140, 66), (145, 66), (150, 66)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 1), (85, 1), (90, 2), (95, 2), (100, 5), (105, 11), (110, 11), (115, 15), (120, 19), (125, 19), (130, 19), (135, 15), (140, 14), (145, 14), (150, 14)]。division_mode=stable_no_division, daughter_size_ratio=None。R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs 多断片化）を特定するための決定的対照。
