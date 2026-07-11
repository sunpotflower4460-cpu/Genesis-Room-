# NOTES — hyp-task3-3-high-Dv

## 出たもの（measured_by / detected）
```json
{
  "detected": {},
  "measured_by": {
    "prediction": "D_vを上げると核形成なし窓が成長域(低比)へ広がる（H-B1安全比が下がる）",
    "verdict": "falsified",
    "points": [
      {
        "kp": 0.03,
        "kd": 0.04,
        "D_v": 2.0,
        "nucleated": true,
        "note": "H-B1既知(再利用)"
      },
      {
        "kp": 0.03,
        "kd": 0.04,
        "D_v": 6.0,
        "nucleated": true,
        "v_background_max": 0.3,
        "u_local_max": 0.935
      },
      {
        "kp": 0.03,
        "kd": 0.04,
        "D_v": 12.0,
        "nucleated": true,
        "v_background_max": 0.3,
        "u_local_max": 0.935
      },
      {
        "kp": 0.03,
        "kd": 0.05,
        "D_v": 2.0,
        "nucleated": true,
        "note": "H-B1既知(再利用)"
      },
      {
        "kp": 0.03,
        "kd": 0.05,
        "D_v": 6.0,
        "nucleated": true,
        "v_background_max": 0.3,
        "u_local_max": 0.558
      }
    ],
    "mechanistic_note": "D_vはvの拡散(再分配)を速めるが、局所反応速度kp*v(uへの変換速度)自体は変えないため、供給が速くても局所生成が追いつき核形成する。"
  }
}
```

## 出なかったもの・未達
reached_level=None, candidate_level=None

## EXPLORATION 所見
課題3-3(高D_v)。判定=falsified。D_vはvの再分配を速めるが局所反応速度は変えないため無効。
