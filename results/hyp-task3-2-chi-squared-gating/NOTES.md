# NOTES — hyp-task3-2-chi-squared-gating

## 出たもの（measured_by / detected）
```json
{
  "detected": {},
  "measured_by": {
    "prediction": "内部分解と表面生産が釣り合うR_cで液滴が安定し、その直上で清潔な分裂",
    "verdict": "partially_supported (rescues droplet from death, but safe-ratio boundary shifts -- needs recalibration)",
    "modification": "R = kp*v*(1-chi(u)) - kd*chi(u)^2  (元は kd*chi(u); 内部のみへの分解ゲート強化)",
    "R0": 15,
    "kp": 0.03,
    "kd": 0.06,
    "rescued_from_death_at_safe_kd": true,
    "still_eventually_nucleates": true,
    "time_series": [
      {
        "t": 0.0,
        "n": 1,
        "top_sizes": [
          15230
        ]
      },
      {
        "t": 25.0,
        "n": 2,
        "top_sizes": [
          6978,
          437
        ]
      },
      {
        "t": 50.0,
        "n": 3,
        "top_sizes": [
          13628,
          6120,
          251
        ]
      },
      {
        "t": 75.0,
        "n": 54,
        "top_sizes": [
          18206,
          11525,
          5005
        ]
      },
      {
        "t": 100.0,
        "n": 63,
        "top_sizes": [
          15439,
          9487,
          4388
        ]
      },
      {
        "t": 125.0,
        "n": 74,
        "top_sizes": [
          14828,
          8611,
          4105
        ]
      }
    ],
    "mass_drift": 0.0,
    "interpretation": "chi^2ゲートは全体的な分解強度を弱めるため、kd=0.06(線形chiでは全R0が死ぬ)でもR0=15が生存・成長した＝「成長を殺さず安全域を使う」方向への前進。ただし背景も結局核形成する(t~75で急増)ため、chi^2モデルにはH-B1と同様の安全window再較正が必要（未実施、時間制約による今回の限界）。"
  }
}
```

## 出なかったもの・未達
reached_level=None, candidate_level=None

## EXPLORATION 所見
課題3-2(Zwicker厳密バランス、chi^2ゲート)。判定=partially_supported (rescues droplet from death, but safe-ratio boundary shifts -- needs recalibration)。有望な部分的成功だが、時間制約で新しい安全windowの再較正(H-B1相当のkp/kdスイープ)は未実施——今後の重要な次の一手。
