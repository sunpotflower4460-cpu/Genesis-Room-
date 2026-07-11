# NOTES — dividing-protocell-3d-probe-R0-09-seed0001

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
    "mass_drift": -2e-15,
    "R0": 9,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 19,
    "division_time": 95.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 40,
    "significant_max_count": 40,
    "raw_final_count": 92,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      2978.0,
      2035.0,
      2033.0,
      2031.0,
      2027.0,
      2027.0,
      2022.0,
      594.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=9.0 probe, t_final=150, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 1), (85, 1), (90, 16), (95, 8), (100, 8), (105, 36), (110, 28), (115, 76), (120, 76), (125, 74), (130, 72), (135, 70), (140, 82), (145, 89), (150, 92)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 1), (85, 1), (90, 1), (95, 5), (100, 8), (105, 9), (110, 17), (115, 20), (120, 20), (125, 20), (130, 32), (135, 32), (140, 32), (145, 40), (150, 40)]。division_mode=stable_no_division, daughter_size_ratio=None。R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs 多断片化）を特定するための決定的対照。
