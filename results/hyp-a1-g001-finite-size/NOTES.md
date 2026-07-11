# NOTES — hyp-a1-g001-finite-size

## 出たもの（measured_by / detected）
```json
{
  "detected": {},
  "measured_by": {
    "prediction": "強BF点(1.5,-1.5)はLとともにrho(L)が非ゼロへ収束(有限サイズ効果)。対照(0,0)は全Lでrho->0(真の秩序化)。",
    "verdict": "falsified",
    "L_values": [
      48,
      96,
      128
    ],
    "strong_BF_point": {
      "b": 1.5,
      "c": -1.5,
      "results": [
        {
          "L": 48,
          "b": 1.5,
          "c": -1.5,
          "steps": 7585,
          "dt_used": 0.019775280898876407,
          "diverged": false,
          "tail_defect_counts": [
            0,
            0,
            0,
            0,
            0,
            0,
            0
          ],
          "tail_density_mean": 0.0,
          "tail_density_final": 0.0,
          "final_defect_count": 0
        },
        {
          "L": 96,
          "b": 1.5,
          "c": -1.5,
          "steps": 7585,
          "dt_used": 0.019775280898876407,
          "diverged": false,
          "tail_defect_counts": [
            0,
            0,
            0,
            0,
            0,
            0,
            0
          ],
          "tail_density_mean": 0.0,
          "tail_density_final": 0.0,
          "final_defect_count": 0
        },
        {
          "L": 128,
          "b": 1.5,
          "c": -1.5,
          "steps": 7585,
          "dt_used": 0.019775280898876407,
          "diverged": false,
          "tail_defect_counts": [
            0,
            0,
            0,
            0,
            0,
            0,
            0
          ],
          "tail_density_mean": 0.0,
          "tail_density_final": 0.0,
          "final_defect_count": 0
        }
      ]
    },
    "control_point": {
      "b": 0.0,
      "c": 0.0,
      "results": [
        {
          "L": 48,
          "b": 0.0,
          "c": 0.0,
          "steps": 3000,
          "dt_used": 0.05,
          "diverged": false,
          "tail_defect_counts": [
            62,
            50,
            24,
            0,
            0,
            0,
            0
          ],
          "tail_density_mean": 0.00017567791005291002,
          "tail_density_final": 0.0,
          "final_defect_count": 0
        },
        {
          "L": 96,
          "b": 0.0,
          "c": 0.0,
          "steps": 3000,
          "dt_used": 0.05,
          "diverged": false,
          "tail_defect_counts": [
            1509,
            1447,
            1390,
            1340,
            1294,
            1248,
            1204
          ],
          "tail_density_mean": 0.0015229724702380952,
          "tail_density_final": 0.0013608579282407408,
          "final_defect_count": 1204
        },
        {
          "L": 128,
          "b": 0.0,
          "c": 0.0,
          "steps": 3000,
          "dt_used": 0.05,
          "diverged": false,
          "tail_defect_counts": [
            3653,
            3433,
            3172,
            2976,
            2751,
            2603,
            2493
          ],
          "tail_density_mean": 0.0014360291617257254,
          "tail_density_final": 0.0011887550354003906,
          "final_defect_count": 2493
        }
      ]
    },
    "bf_density_plateau_detected": false,
    "control_orders_at_all_L": false
  }
}
```

## 出なかったもの・未達
reached_level=None, candidate_level=None

## EXPLORATION 所見
H-A1検証。優先順位2位。判定=falsified。強BF点のrho(L): [0.0, 0.0, 0.0]。対照(0,0)のrho(L): [0.00017567791005291002, 0.0015229724702380952, 0.0014360291617257254]。
