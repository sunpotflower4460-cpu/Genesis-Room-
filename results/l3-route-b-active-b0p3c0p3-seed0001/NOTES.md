# NOTES — l3-route-b-active-b0p3c0p3-seed0001

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
    "params": {
      "b": 0.3,
      "c": 0.3
    },
    "t_final": 400.0,
    "grid": [
      48,
      48,
      48
    ],
    "defect_dynamics": {
      "tail_mean": 0.0,
      "still_coarsening": false,
      "active": false,
      "extinct": true,
      "trend_relative": 0.0,
      "residual_relative_std": 0.0
    },
    "early_phase_mean_coherence": 0.0799,
    "late_phase_mean_coherence": 0.0,
    "coherence_threshold": 0.5,
    "late_phase_judged_coherent": false,
    "coherent_l3_gate": {
      "com_velocity_nonzero": true,
      "circulation_nonzero": false,
      "motion_active_not_extinct": false,
      "structure_coherent_while_moving": false,
      "all_four_pass": false
    },
    "coherence_series_sampled": [
      {
        "t": 6.7,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 26.6,
        "shift": [
          0,
          0,
          0
        ],
        "peak_coherence": 0.0791,
        "n_defects": 619
      },
      {
        "t": 46.6,
        "shift": [
          0,
          0,
          0
        ],
        "peak_coherence": 0.1867,
        "n_defects": 441
      },
      {
        "t": 66.5,
        "shift": [
          1,
          0,
          -1
        ],
        "peak_coherence": 0.1337,
        "n_defects": 364
      },
      {
        "t": 86.5,
        "shift": [
          1,
          0,
          -1
        ],
        "peak_coherence": 0.0811,
        "n_defects": 245
      },
      {
        "t": 106.4,
        "shift": [
          1,
          0,
          -2
        ],
        "peak_coherence": 0.1735,
        "n_defects": 100
      },
      {
        "t": 126.4,
        "shift": [
          0,
          0,
          2
        ],
        "peak_coherence": 0.1505,
        "n_defects": 46
      },
      {
        "t": 146.3,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 166.2,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 186.2,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 206.2,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 226.1,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 246.1,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 266.0,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 285.9,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 305.9,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 325.9,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 345.8,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 365.8,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 385.7,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      },
      {
        "t": 400.0,
        "shift": null,
        "peak_coherence": 0.0,
        "n_defects": 0
      }
    ],
    "defect_count_final": 0,
    "honest_note": "coherenceは位相相関ピーク(隣接スナップショット間の最良剛体シフトとの一致度)で測る代理指標であり、個々の渦芯へのID付与による厳密なobject-trackingではない。ピークが低い/広い場合、パターンが剛体並進でなく再編成(churn)したことを意味する。"
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
L3-ルートB [動くregime候補: 短時間ではactive=True(defect_tail_mean=303.9)。長時間でのcoherence(構造保持)を実測する。] (b,c)=(0.30,0.30), t_final=400。defect: 48683->0 (active=False, still_coarsening=False, extinct=True)。coherence: early=0.07994333333333332 late=0.0 (threshold=0.50, late_coherent=False)。coherent_l3_gate全通過=False。動くregime候補: 短時間ではactive=True(defect_tail_mean=303.9)。長時間でのcoherence(構造保持)を実測する。
