# NOTES — dividing-protocell-3d-probe-R0-13-seed0001

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
    "mass_drift": 2e-15,
    "R0": 13,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 16,
    "division_time": 80.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 25,
    "significant_max_count": 25,
    "raw_final_count": 117,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      1349.0,
      1344.0,
      1342.0,
      1339.0,
      1332.0,
      1331.0,
      1284.0,
      1175.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=13.0 probe, t_final=150, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 2), (20, 2), (25, 2), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 8), (85, 8), (90, 11), (95, 25), (100, 50), (105, 30), (110, 51), (115, 51), (120, 59), (125, 62), (130, 62), (135, 68), (140, 106), (145, 106), (150, 117)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 2), (85, 2), (90, 2), (95, 2), (100, 5), (105, 6), (110, 6), (115, 6), (120, 9), (125, 9), (130, 9), (135, 9), (140, 20), (145, 20), (150, 25)]。division_mode=stable_no_division, daughter_size_ratio=None。R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs 多断片化）を特定するための決定的対照。
