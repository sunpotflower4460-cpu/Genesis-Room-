# NOTES — hyp-b2-phase2-fine-sweep

## 出たもの（measured_by / detected）
```json
{
  "detected": {},
  "measured_by": {
    "prediction": "R0=14-17付近でより対称的な分裂(比->0.5)、方針転換理由はモジュールdocstring参照",
    "pivot_reason": "H-B1のkd/kp安全window(比>=2.0)はPhase1で全R0域(4-24)にわたり液滴を消滅させる regime と判明。growth-permittingなkd=0.04に戻し、既存の突発バースト検出器で交絡排除する方針に転換。",
    "verdict": "supported",
    "phase": "2 (fine sweep near R_c)",
    "R0_sweep": [
      {
        "R0": 12,
        "mode": "multi_fragmentation",
        "genuine_division": true,
        "daughter_size_ratio": 0.4273,
        "division_time": 85.0,
        "significant_final_count": 25,
        "significant_max_count": 25,
        "raw_final_count": 114,
        "burst_detected": false,
        "mass_drift": -2.220446049250313e-16,
        "diverged": false
      },
      {
        "R0": 13,
        "mode": "stable_no_division",
        "genuine_division": false,
        "daughter_size_ratio": null,
        "division_time": 80.0,
        "significant_final_count": 25,
        "significant_max_count": 25,
        "raw_final_count": 117,
        "burst_detected": true,
        "mass_drift": 2.3869795029440866e-15,
        "diverged": false
      },
      {
        "R0": 14,
        "mode": "stable_no_division",
        "genuine_division": false,
        "daughter_size_ratio": null,
        "division_time": 80.0,
        "significant_final_count": 14,
        "significant_max_count": 14,
        "raw_final_count": 68,
        "burst_detected": true,
        "mass_drift": -1.7763568394002505e-15,
        "diverged": false
      },
      {
        "R0": 15,
        "mode": "multi_fragmentation",
        "genuine_division": true,
        "daughter_size_ratio": 0.4397,
        "division_time": 75.0,
        "significant_final_count": 20,
        "significant_max_count": 20,
        "raw_final_count": 109,
        "burst_detected": false,
        "mass_drift": 6.661338147750939e-16,
        "diverged": false
      },
      {
        "R0": 16,
        "mode": "multi_fragmentation",
        "genuine_division": true,
        "daughter_size_ratio": 0.672,
        "division_time": 75.0,
        "significant_final_count": 11,
        "significant_max_count": 11,
        "raw_final_count": 63,
        "burst_detected": false,
        "mass_drift": 0.0,
        "diverged": false
      },
      {
        "R0": 17,
        "mode": "stable_no_division",
        "genuine_division": false,
        "daughter_size_ratio": null,
        "division_time": 15.0,
        "significant_final_count": 10,
        "significant_max_count": 11,
        "raw_final_count": 62,
        "burst_detected": true,
        "mass_drift": -6.661338147750939e-16,
        "diverged": false
      },
      {
        "R0": 18,
        "mode": "stable_no_division",
        "genuine_division": false,
        "daughter_size_ratio": null,
        "division_time": 15.0,
        "significant_final_count": 10,
        "significant_max_count": 11,
        "raw_final_count": 66,
        "burst_detected": true,
        "mass_drift": -2.220446049250313e-16,
        "diverged": false
      },
      {
        "R0": 19,
        "mode": "stable_no_division",
        "genuine_division": false,
        "daughter_size_ratio": null,
        "division_time": 15.0,
        "significant_final_count": 10,
        "significant_max_count": 11,
        "raw_final_count": 60,
        "burst_detected": true,
        "mass_drift": -1.4432899320127035e-15,
        "diverged": false
      },
      {
        "R0": 20,
        "mode": "multi_fragmentation",
        "genuine_division": true,
        "daughter_size_ratio": 0.1904,
        "division_time": 15.0,
        "significant_final_count": 10,
        "significant_max_count": 10,
        "raw_final_count": 74,
        "burst_detected": false,
        "mass_drift": 1.1102230246251565e-15,
        "diverged": false
      }
    ],
    "best_division": {
      "R0": 16,
      "mode": "multi_fragmentation",
      "genuine_division": true,
      "daughter_size_ratio": 0.672,
      "division_time": 75.0,
      "significant_final_count": 11,
      "significant_max_count": 11,
      "raw_final_count": 63,
      "burst_detected": false,
      "mass_drift": 0.0,
      "diverged": false
    }
  }
}
```

## 出なかったもの・未達
reached_level=None, candidate_level=None

## EXPLORATION 所見
H-B2 Phase2。方針転換: H-B1安全windowが成長を殺すため growth-permitting kd=0.04 + 既存burst検出器を使用（理由はdocstring）。判定=supported。
