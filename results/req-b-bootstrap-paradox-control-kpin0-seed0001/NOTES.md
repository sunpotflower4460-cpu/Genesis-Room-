# NOTES — req-b-bootstrap-paradox-control-kpin0-seed0001

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
    "self_maintaining_closure": true,
    "growth_division_inheritance": false,
    "selection_open_ended": false
  },
  "measured_by": {
    "u_bg_initial": 0.12029067695966283,
    "u_bg_peak": 0.4115177126277038,
    "u_bg_peak_time": 91.4,
    "u_bg_final": 0.1776709969592082,
    "u_sp_low": 0.21132486540518713,
    "hump_shape_confirmed": true,
    "background_clean_at_end": true,
    "originated": true,
    "final_n_droplets": 181,
    "mass_drift": -4.625189120588402e-13,
    "u_bg_series_sampled": [
      0.12029067695966283,
      0.1660072227399366,
      0.20792333749291528,
      0.24616185705542287,
      0.28074606723888196,
      0.3115749805062212,
      0.33842761823272954,
      0.36104093169955787,
      0.37927211565433555,
      0.393260476792453,
      0.40343913998388276,
      0.4102157262639526,
      0.4059116493744918,
      0.3572783948320233,
      0.27186553815681685,
      0.2148299551992021,
      0.19422944844143974,
      0.18662822403073762,
      0.1833270763894912,
      0.18168520614511344,
      0.18072293663938022,
      0.1801437785000909,
      0.17953671820548667,
      0.17898344151711162,
      0.1783941670754232
    ]
  }
}
```

## 出なかったもの・未達
reached_level=6, candidate_level=7

## EXPLORATION 所見
依頼B [検証] u_mean=0.12, t_final=200, grid=(96, 96, 96), seed=1。u_bg山形: initial=0.1203 peak=0.4115(t=91.4) final=0.1777 (sp_low=0.2113)。決定的対照: k_p_in=0（区画化なし、種のみ）。背景が汚れたまま（掃除されない）ことを期待。
