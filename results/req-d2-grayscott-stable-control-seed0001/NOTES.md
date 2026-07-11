# NOTES — req-d2-grayscott-stable-control-seed0001

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
    "label": "stable_region_control",
    "F": 0.022,
    "k": 0.051,
    "n_spots_initial": 1,
    "n_spots_max": 8,
    "n_spots_final": 2,
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
        2
      ],
      [
        3200.0,
        1
      ],
      [
        3400.0,
        6
      ],
      [
        3600.0,
        5
      ],
      [
        3800.0,
        1
      ],
      [
        4000.0,
        1
      ],
      [
        4200.0,
        8
      ],
      [
        4400.0,
        1
      ],
      [
        4600.0,
        8
      ],
      [
        4800.0,
        1
      ],
      [
        5000.0,
        1
      ],
      [
        5200.0,
        4
      ],
      [
        5400.0,
        3
      ],
      [
        5600.0,
        2
      ],
      [
        5800.0,
        3
      ],
      [
        5999.0,
        2
      ]
    ],
    "sustained_growth": false,
    "bounded_oscillation": true,
    "honest_caveat": "n_spots系列は単調増加ではなく1と最大8の間を反復振動し(t=2800までの時間の約47パーセントはn=1のまま休眠)、最終的にn=2に落ち着く(最大値の25パーセントまで低下)。mitosis域(F=0.0367,k=0.0649)は単調に増加しn=30-35で安定するのと対照的。単純なreplicated判定(n_max大なる1かつn_final>=2)は両者を同じに分類するが、生データは質的に異なる: mitosisは持続的な自己複製、こちらは境界のある一時的な分岐・再融合(sustained_growthでない)。falsification条件の「対照でも同様に分裂する」を字義通りには満たすが、その中身(持続的増殖 vs 有界振動)は明確に異なるため、パラメータ依存性という主張は生データレベルでは支持される(部分訂正)。"
  }
}
```

## 出なかったもの・未達
reached_level=7, candidate_level=8

## EXPLORATION 所見
依頼D [訂正/精査] 生データ精査により訂正を追記。n_spots系列は単調増加ではなく1と最大8の間を反復振動し(t=2800までの時間の約47パーセントはn=1のまま休眠)、最終的にn=2に落ち着く(最大値の25パーセントまで低下)。mitosis域(F=0.0367,k=0.0649)は単調に増加しn=30-35で安定するのと対照的。単純なreplicated判定(n_max大なる1かつn_final>=2)は両者を同じに分類するが、生データは質的に異なる: mitosisは持続的な自己複製、こちらは境界のある一時的な分岐・再融合(sustained_growthでない)。falsification条件の「対照でも同様に分裂する」を字義通りには満たすが、その中身(持続的増殖 vs 有界振動)は明確に異なるため、パラメータ依存性という主張は生データレベルでは支持される(部分訂正)。
