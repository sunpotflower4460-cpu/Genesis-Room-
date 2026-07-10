# NOTES — g003-model-h-3d-seed0001

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
    "variance_growth": 6e-06,
    "structure_factor_peak_k": 3.0,
    "structure_factor_prominence": 14.773831,
    "correlation_length": 0.333333,
    "defect_count": 0,
    "kinetic_energy": 0.0,
    "circulation_proxy": 2.3e-05,
    "coarsening_length_initial": 0.40855960083028076,
    "coarsening_length_final": 1.7640279751769552,
    "coarsening_exponent_proxy": -0.001450657680612016,
    "max_kinetic_energy": 1.0498619168060757e-09,
    "final_kinetic_energy": 1.0498619168060757e-09,
    "mass_drift": 8.334089847275919e-12
  }
}
```

## 出なかったもの・未達
reached_level=5, candidate_level=6

## EXPLORATION 所見
[検証] seed 1/2/3 は全て Level 5。C=0 対照では流れが消え、capillary coupling が Lv3/Lv5 の決定因子であることを確認した。

[探索] critical mean=0 は bicontinuous coarsening。off-critical phi_mean=0.25 は component count/length が変わる droplet-like regime 候補（Level 0, L_final=1.965）。高C低nu push は max KE=7.368e-09, exponent proxy=-0.001 で、粗い 32^3 では持続循環セルや Lv5 超えの秩序は未確認。Rayleigh-Plateau/トポロジー遷移は可視的断面だけでは主張せず、数値測定器未実装の frontier として残す。

制限: 依頼書の 128^3 正式格子ではなく監査用 coarse global 3D (24^3-40^3, primary 32^3) で実行。正式昇格前に 96^3/128^3 長時間で再実行が必要。無断 clip なし。
