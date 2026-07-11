# NOTES — dividing-protocell-3d-probe-R0-13-res48-seed0001

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
    "mass_drift": -2.38e-13,
    "R0": 13,
    "division_mode": "stable_no_division",
    "division_snapshot_index": 16,
    "division_time": 80.0,
    "daughter_size_ratio": null,
    "significant_initial_count": 1,
    "significant_final_count": 24,
    "significant_max_count": 24,
    "raw_final_count": 113,
    "background_contamination_detected": true,
    "final_sizes_physical_volume": [
      2808.8888888888882,
      1455.4074074074072,
      1410.37037037037,
      1384.2962962962958,
      1372.4444444444441,
      1360.5925925925922,
      1194.6666666666663,
      1192.296296296296
    ]
  }
}
```

## 出なかったもの・未達
reached_level=4, candidate_level=5

## EXPLORATION 所見
依頼3 [決定的対照/診断、自然発生の主張ではない] R0=13.0 probe, t_final=150, grid=(48, 48, 48), dx=1.333。raw count time series: [(0, 1), (5, 1), (10, 1), (15, 2), (20, 2), (25, 2), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 2), (85, 8), (90, 11), (95, 11), (100, 84), (105, 32), (110, 42), (115, 54), (120, 58), (125, 68), (130, 68), (135, 68), (140, 105), (145, 107), (150, 113)]。significant(>=5%初期体積) count time series: [(0, 1), (5, 1), (10, 1), (15, 1), (20, 1), (25, 1), (30, 1), (35, 1), (40, 1), (45, 1), (50, 1), (55, 1), (60, 1), (65, 1), (70, 1), (75, 1), (80, 2), (85, 2), (90, 2), (95, 2), (100, 5), (105, 6), (110, 6), (115, 8), (120, 9), (125, 9), (130, 9), (135, 9), (140, 20), (145, 21), (150, 24)]。division_mode=stable_no_division, daughter_size_ratio=None。解像度収束チェック: 物理領域 L=64 を固定し N=48 vs N=64(dx=1.0, 上の R0=13 probe) で比較。
