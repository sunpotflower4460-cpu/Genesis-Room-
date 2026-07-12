# NOTES — l3-route-a-subcritical-Ra300-seed0001

## 出たもの（measured_by / detected）
```json
{
  "detected": {
    "difference": true,
    "localization": false,
    "spontaneous_motion": false,
    "circulation": false,
    "persistent_individuality": false,
    "co_differentiation": false,
    "self_maintaining_closure": false,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "Ra": 300.0,
    "Ra_c_theory": 657.5113644795163,
    "Pr": 1.0,
    "grid": [
      17,
      24,
      24
    ],
    "t_final": 1.5,
    "dt": 0.0005,
    "diverged": false,
    "w_rms_series": [
      2.3e-05,
      9.3e-05,
      6.3e-05,
      4.4e-05,
      3.1e-05,
      2.3e-05,
      1.7e-05,
      1.2e-05,
      9e-06,
      7e-06,
      5e-06,
      4e-06,
      3e-06,
      2e-06,
      2e-06,
      1e-06,
      1e-06,
      1e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "circulation_proxy_series": [
      4e-06,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0,
      0.0
    ],
    "fitted_growth_rate": -4.583302057298952,
    "growing": false,
    "circulation_nonzero": false,
    "coherence_series": [
      0.3831,
      0.2859,
      0.1745,
      0.1195,
      0.0976,
      0.0741,
      0.0696,
      0.0605,
      0.0568,
      0.0505,
      0.043,
      0.0376,
      0.0351,
      0.0348,
      0.0347,
      0.0347,
      0.0347,
      0.0331,
      0.0278,
      0.023,
      0.0208,
      0.0187,
      0.0166,
      0.0149
    ],
    "mean_coherence": 0.0742,
    "early_phase_mean_coherence": 0.1211,
    "late_phase_mean_coherence": 0.0274,
    "coherent_gate_pass": false,
    "l3_gate": {
      "com_velocity_nonzero": false,
      "circulation_nonzero": false,
      "coherent_while_growing": false
    },
    "honest_limitation": "axis順序バグ修正後、Ra=1200・t_final=3で数値発散は解消し、非線形飽和(w_rmsの伸びの明確な頭打ち)まで到達を確認した。ただし完全な定常状態(w_rmsが時間的に一定)への収束は本run終了時点では確認できておらず(循環proxyはまだ緩やかに増加中)、『飽和しつつある』段階の観測である。coherenceは位相相関の全区間平均でなく後期(構造形成後)半分の平均で判定する——初期ノイズ遷移は多モード競合でcoherenceが見かけ上低くなるため(ルートBと同じ手法)。"
  }
}
```

## 出なかったもの・未達
reached_level=2, candidate_level=3

## EXPLORATION 所見
L3-ルートA [決定的対照(亜臨界): 減衰、動かないことを期待。] Ra=300 (Ra_c理論値=657.5), t_final=1.50, grid=(17, 24, 24)。w_rms: ['2.28e-05', '9.32e-05', '6.27e-05', '4.38e-05', '3.13e-05', '2.27e-05', '1.66e-05', '1.22e-05', '9.03e-06', '6.69e-06', '4.97e-06', '3.70e-06', '2.75e-06', '2.05e-06', '1.53e-06', '1.14e-06', '8.50e-07', '6.34e-07', '4.72e-07', '3.52e-07', '2.63e-07', '1.96e-07', '1.46e-07', '1.09e-07', '8.17e-08']。fitted_growth_rate=-4.583302057298952, growing=False, circulation_nonzero=False, mean_coherence=0.07423653987759694 (early=0.12105677772339822, late=0.02741630203179569, coherent_gate=False)。diverged=False。決定的対照(亜臨界): 減衰、動かないことを期待。
