# NOTES — vortex-a-opppair-dissipative-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": true,
    "spontaneous_motion": false,
    "circulation": true,
    "persistent_individuality": false,
    "co_differentiation": false,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "config": "opposite_pair",
    "gamma": 0.3,
    "n_cores_initial": 2,
    "n_cores_final": 0,
    "coexisted": false,
    "annihilated": true,
    "net_winding_series": [
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0,
      0
    ],
    "winding_conserved": true,
    "energy0": -15684.429,
    "energy_final": -16378.498,
    "energy_drift_frac": -0.04425,
    "norm0": 32628.658,
    "norm_final": 32757.383,
    "core_motion": {
      "translated": null,
      "cum_rotation_deg": null
    },
    "n_cores_series": [
      [
        0.0,
        2
      ],
      [
        1.2,
        2
      ],
      [
        2.4,
        2
      ],
      [
        3.6,
        2
      ],
      [
        4.8,
        2
      ],
      [
        6.0,
        2
      ],
      [
        7.2,
        2
      ],
      [
        8.4,
        2
      ],
      [
        9.6,
        0
      ],
      [
        10.8,
        0
      ],
      [
        12.0,
        0
      ],
      [
        13.2,
        0
      ],
      [
        14.4,
        0
      ],
      [
        15.6,
        0
      ],
      [
        16.8,
        0
      ],
      [
        18.0,
        0
      ],
      [
        19.2,
        0
      ],
      [
        20.4,
        0
      ],
      [
        21.6,
        0
      ],
      [
        22.8,
        0
      ],
      [
        23.98,
        0
      ]
    ],
    "min_density_final": 0.9884
  }
}
```

## 出なかったもの・未達
reached_level=3, candidate_level=4

## EXPLORATION 所見
渦-依頼A [opposite_pair] gamma=0.30, t_final=24, grid=(64, 64, 64)。n_cores 2->0, coexisted=False, annihilated=True, winding_conserved=True, E_drift=-0.0443, motion={'translated': None, 'cum_rotation_deg': None}。決定的対照: 同一配置でgammaのみ0→0.3。対消滅の予言。
