# NOTES — dividing-protocell-3d-probe-R0-17-longrun-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": false,
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
    "mass_drift": 1e-15,
    "R0": 17,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 8,
    "division_time": 80.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 0,
    "significant_max_count": 11,
    "raw_final_count": 139,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      1078.0,
      1075.0,
      1074.0,
      1070.0,
      1066.0,
      1066.0,
      911.0,
      620.0
    ]
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=17.0 probe, t_final=300, grid=(64, 64, 64), dx=1.000。raw count time series: [(0, 1), (10, 1), (20, 2), (30, 2), (40, 2), (50, 2), (60, 2), (70, 21), (80, 3), (90, 38), (100, 45), (110, 63), (120, 47), (130, 47), (140, 58), (150, 62), (160, 79), (170, 100), (180, 121), (190, 121), (200, 133), (210, 136), (220, 136), (230, 136), (240, 136), (250, 137), (260, 139), (270, 139), (280, 139), (290, 139), (300, 139)]。significant(>=5%初期体積) count time series: [(0, 1), (10, 1), (20, 1), (30, 1), (40, 1), (50, 1), (60, 1), (70, 1), (80, 2), (90, 8), (100, 11), (110, 11), (120, 11), (130, 11), (140, 10), (150, 10), (160, 9), (170, 9), (180, 6), (190, 6), (200, 6), (210, 6), (220, 6), (230, 6), (240, 6), (250, 6), (260, 6), (270, 6), (280, 6), (290, 1), (300, 0)]。division_mode=stable_no_division, daughter_size_ratio=None。[探索] grow-divide サイクル探索: t_final を2倍に延長し、分裂後の娘液滴が再成長・再分裂するかを追跡。
