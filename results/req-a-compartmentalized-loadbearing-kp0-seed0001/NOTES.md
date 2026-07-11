# NOTES — req-a-compartmentalized-loadbearing-kp0-seed0001

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
    "u_mean": 0.35,
    "variance_growth": 0.003453,
    "structure_factor_prominence": 17.927,
    "final_n_droplets": 110,
    "final_top_sizes": [
      128.0,
      123.0,
      114.0,
      106.0,
      104.0
    ],
    "u_bg_final": 0.1556393793903142,
    "u_bg_max": 0.34996913348372155,
    "u_sp_low": 0.21132486540518713,
    "background_clean": true,
    "mass_drift": -4.4941828036826337e-13,
    "load_bearing_correction": "自動判定はn_droplets>=1のみでpersistent_individualityをTrueとしていたが、これは相分離(spinodal decomposition、反応と無関係にCH力学だけで起きる)による核形成と、反応駆動の大きな支配的個体形成を区別できていなかった。正しい比較: 主実験(kp=0.05)の最大液滴=273753 voxel（ドメインの31%）に対し、この対照(kp=0)は最大液滴=128 voxel（比=0.0005、約2139xの差）で、110個の小さな静止/縮小中の液滴に留まった——支配的な個体は形成されなかった。2D確認値(total_u 4744->0で死亡)と定性的に整合する、正しいload-bearing死の確認。生産(kp)は、液滴が核形成すること自体には不要（相分離が担う）が、核形成した液滴が大きな自己維持個体へconsolidateするには本質的——これが正確な結論。"
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
# NOTES — req-a-compartmentalized-loadbearing-kp0-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": true,
    "spontaneous_motion": false,
    "circulation": false,
    "persistent_individuality": true,
    "co_differentiation": false,
    "self_maintaining_closure": true,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "u_mean": 0.35,
    "variance_growth": 0.003453,
    "structure_factor_prominence": 17.927,
    "final_n_droplets": 110,
    "final_top_sizes": [
      128.0,
      123.0,
      114.0,
      106.0,
      104.0
    ],
    "u_bg_final": 0.1556393793903142,
    "u_bg_max": 0.34996913348372155,
    "u_sp_low": 0.21132486540518713,
    "background_clean": true,
    "mass_drift": -4.4941828036826337e-13
  }
}
```

## 出なかったもの・未達
reached_level=6, candidate_level=7

## EXPLORATION 所見
依頼A [検証] u_mean=0.35, t_final=150, grid=(96, 96, 96), seed=1。0≠無 climb: reached_level=6。u_bg: [0.34996913348372155, 0.3458950171842863, 0.3421146681857205, 0.3378238931988989, 0.30686084555310233, 0.2391580617207677, 0.19822650822311036, 0.17959757529719192, 0.16936459297498277, 0.16325081254707322, 0.1594846227140988, 0.15736229894779766, 0.1562590504902828, 0.15574216905129296, 0.15563044197188097]。決定的対照2 (load-bearing): k_p=0 で生産を切る。個体が死ねば代謝が本質=Lv6が本物であることの証拠。


[訂正] load-bearing死の判定を修正。詳細はmeasured_by.load_bearing_correction参照。reached_level=2(相分離のみ、大きな個体化には反応が必須)、role=N（負の結果、求める現象=支配的個体化が起きないことを確認）。
