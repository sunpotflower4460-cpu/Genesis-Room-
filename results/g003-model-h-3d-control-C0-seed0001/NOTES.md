# NOTES — g003-model-h-3d-control-C0-seed0001

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
    "variance_growth": 0.028243,
    "structure_factor_peak_k": 4.0,
    "structure_factor_prominence": 19.401684,
    "correlation_length": 0.25,
    "defect_count": 0,
    "kinetic_energy": 0.0,
    "circulation_proxy": 0.0,
    "occupied_fraction_phase_pos": 0.500868,
    "connected_components_pos": 2,
    "connected_components_neg": 2,
    "coarsening_rate_exponent": 0.1261875960201641,
    "coarsening_length_final": 1.592054,
    "final_kinetic_energy": 0.0
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
依頼1 [検証]: 2D で見つけた Lv1(分散/構造因子)->Lv2(界面/連結成分)->Lv3(KE>0)->Lv5(共分化、C=0対照で流れ消失+粗大化減速により相互依存を確認) が 3D 48^3 grid, t_final=100 (steps=2000, dt=0.050) で立つかを検証。mode=coarse-global-3d 相当（依頼書の 128^3 full-3d には未到達、下記注意点を参照）。
