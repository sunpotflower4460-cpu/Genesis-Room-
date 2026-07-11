# NOTES — g003-model-h-3d-seed0003

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": true,
    "spontaneous_motion": true,
    "circulation": true,
    "persistent_individuality": false,
    "co_differentiation": true,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "variance_growth": 0.036444,
    "structure_factor_peak_k": 2.0,
    "structure_factor_prominence": 23.28779,
    "correlation_length": 0.5,
    "defect_count": 0,
    "kinetic_energy": 0.005605,
    "circulation_proxy": 0.027186,
    "occupied_fraction_phase_pos": 0.500289,
    "connected_components_pos": 1,
    "connected_components_neg": 2,
    "control_C0_final_kinetic_energy": 0.0,
    "control_C0_coarsening_length": 1.592054,
    "primary_coarsening_length": 2.20179,
    "coarsening_rate_exponent": 0.24351595822631186,
    "coarsening_length_final": 2.20179,
    "final_kinetic_energy": 0.005604661
  }
}
```

## 出なかったもの・未達
reached_level=5, candidate_level=6

## EXPLORATION 所見
依頼1 [検証]: 2D で見つけた Lv1(分散/構造因子)->Lv2(界面/連結成分)->Lv3(KE>0)->Lv5(共分化、C=0対照で流れ消失+粗大化減速により相互依存を確認) が 3D 48^3 grid, t_final=100 (steps=2000, dt=0.050) で立つかを検証。mode=coarse-global-3d 相当（依頼書の 128^3 full-3d には未到達、下記注意点を参照）。
