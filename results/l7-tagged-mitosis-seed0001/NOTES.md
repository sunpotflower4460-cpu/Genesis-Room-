# NOTES — l7-tagged-mitosis-seed0001

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
    "use_tags": true,
    "F": 0.0367,
    "k": 0.0649,
    "n_seeds": 7,
    "founder_tags": [
      1,
      1,
      0,
      0,
      0,
      1,
      1
    ],
    "grid": [
      64,
      64,
      64
    ],
    "steps": 6000,
    "n_spots_initial": 8,
    "n_spots_max": 43,
    "n_spots_final": 43,
    "n_series_sampled": [
      [
        0.0,
        8
      ],
      [
        200.0,
        9
      ],
      [
        400.0,
        9
      ],
      [
        600.0,
        10
      ],
      [
        800.0,
        14
      ],
      [
        1000.0,
        15
      ],
      [
        1200.0,
        18
      ],
      [
        1400.0,
        20
      ],
      [
        1600.0,
        22
      ],
      [
        1800.0,
        24
      ],
      [
        2000.0,
        25
      ],
      [
        2200.0,
        28
      ],
      [
        2400.0,
        29
      ],
      [
        2600.0,
        33
      ],
      [
        2800.0,
        38
      ],
      [
        3000.0,
        38
      ],
      [
        3200.0,
        41
      ],
      [
        3400.0,
        39
      ],
      [
        3600.0,
        40
      ],
      [
        3800.0,
        39
      ],
      [
        4000.0,
        38
      ],
      [
        4200.0,
        37
      ],
      [
        4400.0,
        36
      ],
      [
        4600.0,
        38
      ],
      [
        4800.0,
        39
      ],
      [
        5000.0,
        41
      ],
      [
        5200.0,
        41
      ],
      [
        5400.0,
        42
      ],
      [
        5600.0,
        42
      ],
      [
        5800.0,
        43
      ],
      [
        5999.0,
        43
      ]
    ],
    "division_not_seeded": true,
    "final_spot_purities": [
      1.0,
      0.9935,
      1.0,
      1.0,
      1.0,
      0.9859,
      0.9978,
      1.0,
      1.0,
      1.0,
      1.0,
      0.9994,
      1.0,
      1.0,
      1.0,
      0.9958,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      0.9778,
      1.0,
      1.0,
      0.6273,
      0.986,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      0.9974,
      0.9942,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0,
      1.0
    ],
    "clean_fraction": 1.0,
    "state_inherited": true,
    "locality_checks": 800,
    "locality_violations": 0,
    "locality_ok": true,
    "volume_rel_change_typical": 0.002817,
    "volume_rel_change_outliers": 0,
    "volume_continuous": true,
    "accounting_warmup_snapshots_excluded": 10,
    "divisions_in_warmup_excluded_from_accounting": 4,
    "accounting_consistent": true,
    "l7_gate": {
      "division_not_seeded": true,
      "state_inherited": true,
      "accounting_consistent": true,
      "reached_L7": true
    },
    "l7_partial": false,
    "honest_note": "state_inheritedは「各スポットのタグがclean(0/1にロック)されているか」の割合(bistable purity>0.6の閾値)で測る——依頼書の定義通り。娘スポットが親と同じ系統IDを持つことの明示的なlineage trackingは行っていない(親子関係の直接追跡でなく、事後のタグ純度で代理測定)。accounting_consistentは(a)新規スポットが直前の既存スポット近傍(周期距離)にのみ出現するか(空の場所からの湧き出しでないか)、(b)スポット総体積の相対変化に分裂時の異常な飛躍がないか、の2条件。(b)の「通常時」基準は創始者スポットが安定サイズへ整定するまでの初期過渡(warmup、accounting_warmup_snapshots_excluded)を除外して求める——この間の体積変化は自己複製でなく初期条件からの物理的な整定であり、除外しないと真の分裂ジャンプの基準がゆがむ。同じ理由でwarmup中に起きた分裂(divisions_in_warmup_excluded_from_accounting)はaccounting評価の対象外とする(整定と複製が混ざり評価不能なため)。"
  }
}
```

## 出なかったもの・未達
reached_level=7, candidate_level=8

## EXPLORATION 所見
確認B [タグあり: 分裂+遺伝+会計の全3要件を期待。] F=0.0367 k=0.0649, steps=6000, n_seeds=7, founder_tags=[1, 1, 0, 0, 0, 1, 1]。n_spots: 8->43(max=43)。division_not_seeded=True, clean_fraction=1.0, state_inherited=True, locality_ok=True(0/800), volume_continuous=True(outliers=0), accounting_consistent=True -> reached_L7=True(l7_partial=False)。タグあり: 分裂+遺伝+会計の全3要件を期待。
