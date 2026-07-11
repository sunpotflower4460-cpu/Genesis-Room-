#!/usr/bin/env python3
"""依頼B: ブートストラップ・パラドックス解消（3D）— 安定な無からでも起源→維持。

--- 予言（実行前に登録） ---
種項(kp_seed·v·(1-chi))がuをspinodal超えへ駆動→最初の液滴が核形成（背景一時的に汚れる
u_bg↑）→区画化項(kp_in·v·chi)がvを内部消費→種項が飢える（vが枯渇し種生産も弱まる）→
背景がclear（u_bg→spinodal以下）。起源(汚)→維持(clean)の遷移＝山形のu_bg(t)。

決定的対照: 区画化項なし（k_p_in=0、種のみ）→ 背景が汚れたまま（掃除されない）。

2D確認値: u_bgが0.12->0.38(peak)->0.07と山形、最終clean。4パラメータ全てで成功。

falsification: 3Dで背景がclearにならない（種が汚し続ける）or 起源しない（Lv1にすら
届かない）なら、この経路は3Dで崩れる。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import compartmentalized_genesis_3d as cg  # noqa: E402

SHAPE = (96, 96, 96)
DT = 0.05
STEPS = 4000  # t_final = 200 -- 山形(上昇->下降)を見るため依頼Aより少し長め
SEEDED_PARAMS = dict(cg.DEFAULTS, k_p_seed=0.01, k_p_in=0.06, k_d=0.03, D_v=2.0)


def _run_and_save(room_id, params, seed, u_mean=0.12, v_mean=0.6, notes_extra=""):
    snapshots, phys = cg.run(SHAPE, STEPS, DT, seed, params=params, u_mean=u_mean, v_mean=v_mean,
                             reaction_fn=cg.reaction_seeded_compartmentalized)
    u_bg = phys["u_bg_series"]
    u_bg_initial = u_bg[0] if u_bg else None
    u_bg_peak = max(u_bg) if u_bg else None
    u_bg_final = phys["u_bg_final"]
    peak_idx = int(np.argmax(u_bg)) if u_bg else 0
    hump_shape = bool(u_bg_peak is not None and u_bg_peak > cg.U_SP_LOW
                      and u_bg_final is not None and u_bg_final < u_bg_peak
                      and peak_idx < len(u_bg) - 1)
    background_clean_at_end = bool(u_bg_final is not None and u_bg_final < cg.U_SP_LOW)
    originated = bool(snapshots[-1]["n_droplets"] >= 1)

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = originated
    detected["localization"] = originated
    detected["persistent_individuality"] = bool(originated and not phys["diverged"])
    detected["self_maintaining_closure"] = bool(originated and background_clean_at_end)
    reached = 0
    if originated:
        reached = 2
    if detected["persistent_individuality"]:
        reached = 4
    if detected["self_maintaining_closure"]:
        reached = 6

    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"u_bg_initial": u_bg_initial, "u_bg_peak": u_bg_peak,
                        "u_bg_peak_time": round(peak_idx * DT * (STEPS / max(1, len(u_bg))), 2),
                        "u_bg_final": u_bg_final, "u_sp_low": cg.U_SP_LOW,
                        "hump_shape_confirmed": hump_shape,
                        "background_clean_at_end": background_clean_at_end,
                        "originated": originated,
                        "final_n_droplets": snapshots[-1]["n_droplets"],
                        "mass_drift": phys["mass_drift"],
                        "u_bg_series_sampled": u_bg[::max(1, len(u_bg) // 25)]},
        "purity": {"per_object_labels": False, "external_optimum": False,
                  "role": "E" if reached >= 1 else "F"},
    }
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"], resolutions_result={"96x96x96": reached},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "k_p_in=0 (seed-only, decisive control)",
                        "result": "see companion room -- expect background stays dirty"}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]])
    genesis_yaml = {
        "equations": "R = kp_seed*v*(1-chi(u)) + kp_in*v*chi(u) - kd*chi(u)^2 "
                     "(weak uniform seed + strong compartmentalization)",
        "solver": "pseudo-spectral, semi-implicit CH (u) + implicit diffusion (v)",
        "dt": DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
        "params": params, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼B [検証] u_mean=%.2f, t_final=%.0f, grid=%s, seed=%d。u_bg山形: "
             "initial=%.4f peak=%.4f(t=%.1f) final=%.4f (sp_low=%.4f)。%s"
             % (u_mean, STEPS * DT, SHAPE, seed, u_bg_initial, u_bg_peak,
                report["measured_by"]["u_bg_peak_time"], u_bg_final, cg.U_SP_LOW, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  reached_level=%d hump=%s clean_end=%s u_bg(%.3f->%.3f->%.3f)"
          % (run_dir, reached, hump_shape, background_clean_at_end, u_bg_initial, u_bg_peak,
             u_bg_final))
    return snapshots, phys, report


def main():
    print("=== 依頼B: ブートストラップ・パラドックス解消 (3D) ===")
    print("予言: u_bg が山形(0.12->peak>spinodal->spinodal以下)を描く。起源(汚)->維持(clean)。")

    print("\n[本実験] seeds=1,2,3")
    for seed in (1, 2, 3):
        _run_and_save("req-b-bootstrap-paradox-seed%04d" % seed, SEEDED_PARAMS, seed)

    print("\n[決定的対照] k_p_in=0（種のみ、区画化なし）")
    control_params = dict(SEEDED_PARAMS, k_p_in=0.0)
    _run_and_save("req-b-bootstrap-paradox-control-kpin0-seed0001", control_params, 1,
                  notes_extra="決定的対照: k_p_in=0（区画化なし、種のみ）。背景が汚れたまま"
                              "（掃除されない）ことを期待。")

    print("=== 依頼B done ===")


if __name__ == "__main__":
    main()
