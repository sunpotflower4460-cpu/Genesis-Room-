# NOTES — dividing-protocell-3d-probe-R0-17-seed0001

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
    "mass_drift": -1e-15,
    "R0": 17,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 3,
    "division_time": 15.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 10,
    "significant_max_count": 11,
    "raw_final_count": 62,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      1779.0,
      1776.0,
      1776.0,
      1376.0,
      1374.0,
      1373.0,
      1372.0,
      1372.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=17.0 probe, t_final=150, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 2), (20, 2), (25, 2), (30, 2), (35, 2), (40, 2), (45, 2), (50, 2), (55, 2), (60, 2), (65, 2), (70, 21), (75, 9), (80, 3), (85, 21), (90, 38), (95, 39), (100, 45), (105, 63), (110, 63), (115, 47), (120, 47), (125, 47), (130, 47), (135, 47), (140, 58), (145, 60), (150, 62)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 2), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 2), (80, 2), (85, 2), (90, 8), (95, 11), (100, 11), (105, 11), (110, 11), (115, 11), (120, 11), (125, 11), (130, 11), (135, 11), (140, 10), (145, 11), (150, 10)]。division_mode=stable_no_division, daughter_size_ratio=None。R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs 多断片化）を特定するための決定的対照。
