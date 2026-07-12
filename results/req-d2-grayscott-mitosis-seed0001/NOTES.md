# NOTES — req-d2-grayscott-mitosis-seed0001

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
    "self_maintaining_closure": false,
    "growth_division_inheritance": true,
    "selection_open_ended": false
  },
  "measured_by": {
    "label": "mitosis_region",
    "F": 0.0367,
    "k": 0.0649,
    "n_spots_initial": 1,
    "n_spots_max": 35,
    "n_spots_final": 33,
    "replicated_passively": true,
    "n_series_sampled": [
      [
        0.0,
        1
      ],
      [
        200.0,
        1
      ],
      [
        400.0,
        1
      ],
      [
        600.0,
        1
      ],
      [
        800.0,
        1
      ],
      [
        1000.0,
        1
      ],
      [
        1200.0,
        1
      ],
      [
        1400.0,
        1
      ],
      [
        1600.0,
        1
      ],
      [
        1800.0,
        1
      ],
      [
        2000.0,
        1
      ],
      [
        2200.0,
        1
      ],
      [
        2400.0,
        1
      ],
      [
        2600.0,
        1
      ],
      [
        2800.0,
        1
      ],
      [
        3000.0,
        1
      ],
      [
        3200.0,
        1
      ],
      [
        3400.0,
        1
      ],
      [
        3600.0,
        1
      ],
      [
        3800.0,
        2
      ],
      [
        4000.0,
        4
      ],
      [
        4200.0,
        9
      ],
      [
        4400.0,
        12
      ],
      [
        4600.0,
        17
      ],
      [
        4800.0,
        23
      ],
      [
        5000.0,
        29
      ],
      [
        5200.0,
        31
      ],
      [
        5400.0,
        34
      ],
      [
        5600.0,
        35
      ],
      [
        5800.0,
        30
      ],
      [
        5999.0,
        33
      ]
    ]
  }
}
```

## 出なかったもの・未達
reached_level=7, candidate_level=8

## EXPLORATION 所見
依頼D [mitosis_region] F=0.0367 k=0.0649, t_final=6000, grid=(64, 64, 64), seed=1。n_spots: 1->33(max=35)。replicated_passively=True。series: [(0.0, 1), (400.0, 1), (800.0, 1), (1200.0, 1), (1600.0, 1), (2000.0, 1), (2400.0, 1), (2800.0, 1), (3200.0, 1), (3600.0, 1), (4000.0, 4), (4400.0, 12), (4800.0, 23), (5200.0, 31), (5600.0, 35), (5999.0, 33)]。mitosisパラメータ域(既知の自己複製スポット)。
