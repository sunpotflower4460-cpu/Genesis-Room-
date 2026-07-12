#!/usr/bin/env python3
"""L3-ルートA【最有力】: Rayleigh-Bénard対流をfull-3Dで（清潔なcoherent L3）。

--- 予言（実行前に登録） ---
静止流体(u=v=w=0)＋Tへの微小ノイズのみ（対流セル・回転・循環は一切初期条件に置かない——
第8監査L3）から、genesis/boussinesq_rb_3d.py（無次元Boussinesq、自由すべり壁＋固定温度差、
domain-doublingトリックで既存の周期擬スペクトル基盤を再利用）を発展させる。

臨界Ra_c(自由すべり) = 27*pi^4/4 ~ 657.51（既知の解析解）を境に、Ra<Ra_cでは静止(摂動が
指数減衰)、Ra>Ra_cでは自発的にw,循環が成長する、と予言する。自己テスト(genesis/
boussinesq_rb_3d.py の __main__、単一モードでの精密な線形安定性検証)で、この定性的な
遷移(部分臨界=減衰、超臨界=成長)を確認済み。

★清潔なcoherent L3の判定: 単なる「動く」でなく、渦線coarseningのような乱流churnと
区別する。growth期のwパターンが空間的に自己相似に増幅する(coherent固有モード的成長)か、
それとも形が絶えず再編成される(turbulent)かを、位相相関(隣接スナップショット間の相関
ピーク)で測る——ルートBと同じcoherence指標を流用。

honest限界: 開発初期、Leray投影(非圧縮性強制)の軸順序バグにより非線形領域で数値発散が
頻発していたが、genesis/boussinesq_rb_3d.py 側でこのバグを発見・修正した(運動学的に
w=0のはずの水平一様モードでwが非ゼロに成長する、というテストで検出)。修正後はRa=1200・
t_final=3で数値発散なく、w_rmsの成長が明確に頭打ちする(非線形飽和)ことを確認した。
ただし本runの終了時点で循環proxy(渦度z成分rms)はまだ緩やかに増加中であり、完全な
定常状態(時間的に一定)への収束そのものは未確認——「飽和しつつある」段階の観測である、
と正直に報告する。またRa=200/1200での自己テスト時の成長・減衰率は理論値と符号は
一致するが大きさは数倍のずれがあり(genesis/boussinesq_rb_3d.py参照)、この定量的な
精度の検証は今回のスコープ外とする。

決定的対照: Ra=300(<Ra_c、静止のはず) vs Ra=1200(>Ra_c、成長・飽和するはず)。

--- 追加の予言（確認A: 4acf1a46-L3_L7_confirm_3D_GenesisRoom_Sonnet5.pdf、実行前に登録） ---
依頼書指定の全要件ゲート(KE≠0 AND 3D循環≠0(common/diagnostics.circulation) AND
Nu(=1+<w*T>)>1 AND late_fluct(KEの後期相対std)<0.05)で、t_final=3では循環proxyが
まだ緩やかに増加中だったため、t_final=8まで延長すれば真の定常状態(late_fluct<0.05)に
到達し、Ra=1200で4要件全通過(coherent L3確定)、Ra=300では全て不通過のまま(静止)、
と予言する。一部の要件のみ通過ならl3_partialとして明記し、reached_level=3には昇格させない。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import boussinesq_rb_3d as rb  # noqa: E402

NZ, NY, NX = 17, 24, 24
H, LX, LY = 1.0, 4.0, 4.0
RA_C_THEORY = 27.0 * np.pi ** 4 / 4.0


def phase_correlation_shift(density0, density1):
    if density0.std() < 1e-12 or density1.std() < 1e-12:
        return None, 0.0
    F0 = np.fft.fftn(density0 - density0.mean())
    F1 = np.fft.fftn(density1 - density1.mean())
    R = F0 * np.conj(F1)
    denom = np.abs(R)
    denom[denom < 1e-12] = 1e-12
    R = R / denom
    corr = np.real(np.fft.ifftn(R))
    peak_idx = np.unravel_index(np.argmax(corr), corr.shape)
    peak_val = float(corr[peak_idx])
    return peak_idx, peak_val


def run_and_save(room_id, Ra, seed=1, t_final=0.6, notes_extra=""):
    p = {"Pr": 1.0, "Ra": Ra}
    dt = 0.0005
    snapshots, phys = rb.run(NZ, NY, NX, H, LX, LY, t_final, seed, params=p,
                             snapshot_every_frac=1.0 / 24, noise_amplitude=0.001, dt=dt)

    w_rms_series = []
    circ_series = []
    coherence_series = []
    ke_series = []
    circ3d_series = []
    nu_series = []
    prev_w = None
    for s in snapshots:
        w_phys = rb.physical_slice(s["w"], NZ)
        u_phys, v_phys = rb.physical_slice(s["u"], NZ), rb.physical_slice(s["v"], NZ)
        T_phys = rb.physical_slice(s["T"], NZ)
        w_rms_series.append(float(np.sqrt(np.mean(w_phys ** 2))))
        # 循環の代理: 水平面内の速度回転(渦度z成分)のrms、中間高度で評価
        mid = NZ // 2
        dv_dx = np.gradient(v_phys[mid], axis=1)
        du_dy = np.gradient(u_phys[mid], axis=0)
        vort_z = dv_dx - du_dy
        circ_series.append(float(np.sqrt(np.mean(vort_z ** 2))))
        # 運動エネルギー・循環(3D渦度、common/diagnostics.circulationを再利用=他Roomと同じ
        # 恣意性のない測定)・Nu(=1+<w*T>、無次元化により伝導フラックスが厳密に1になる規約)
        ke_series.append(diag.kinetic_energy([u_phys, v_phys, w_phys]))
        circ3d_series.append(diag.circulation([u_phys, v_phys, w_phys]))
        nu_series.append(1.0 + float(np.mean(w_phys * T_phys)))
        if prev_w is not None and w_rms_series[-2] > 1e-9:
            _, peak = phase_correlation_shift(prev_w[mid], w_phys[mid])
            coherence_series.append(peak)
        prev_w = w_phys

    times = [s["step"] * phys["dt_used"] for s in snapshots]
    growing = bool(w_rms_series[-1] > 3.0 * max(w_rms_series[0], 1e-12) and not phys["diverged"])
    # 最後の(未発散な)有効フレームまでの指数成長率をフィット
    valid_w = [w for w in w_rms_series if np.isfinite(w) and w > 0]
    if len(valid_w) >= 3:
        t_valid = np.array(times[:len(valid_w)])
        slope = float(np.polyfit(t_valid, np.log(valid_w), 1)[0])
    else:
        slope = None

    mean_coherence = float(np.mean(coherence_series)) if coherence_series else None
    # ルートBと同じ手法: 早期(ノイズ優勢な多モード遷移)と後期(飽和/単一ロールが支配的)に
    # 分けて評価する。全区間平均は初期ノイズ遷移で希釈されるため、「構造が保たれたまま動くか」
    # を問うなら、構造が実際に存在する後期区間の値がより意味を持つ。
    n_half = len(coherence_series) // 2
    early_coh = coherence_series[:n_half] if n_half else []
    late_coh = coherence_series[n_half:] if n_half else []
    early_mean_coherence = float(np.mean(early_coh)) if early_coh else None
    late_mean_coherence = float(np.mean(late_coh)) if late_coh else None
    circulation_nonzero = bool(circ_series[-1] > 2.0 * max(circ_series[0], 1e-12))
    com_velocity_nonzero = growing  # w自体が成長=鉛直流速が非ゼロ

    # --- 依頼(L3/L7指示書)指定の全要件ゲート ---
    # KE(common/diagnostics.kinetic_energyで統一測定、他Roomと同じ関数=恣意性なし)
    ke_nonzero = bool(ke_series[-1] > 3.0 * max(ke_series[0], 1e-15) and not phys["diverged"])
    # 3D循環(common/diagnostics.circulation, common/emergence.pyが使う関数と同一)
    circ3d_nonzero = bool(circ3d_series[-1] > 3.0 * max(circ3d_series[0], 1e-12))
    # 熱流束: Nu=1+<w*T>(無次元化で伝導フラックスがちょうど1になる規約)。後期平均で評価
    n_late = max(1, len(nu_series) // 4)
    nu_late_mean = float(np.mean(nu_series[-n_late:]))
    flux_nonzero = bool(nu_late_mean > 1.02)  # 静止時のノイズ床(~1.00±noise)を明確に超える
    # coherentフラグ: KEのlate-time相対std(late_fluct)<0.05=定常ロール(乱流でない)。
    # ★依頼書が指定する厳密な定義(位相相関ベースのcoherenceとは別の、独立した指標)。
    ke_late = ke_series[-n_late:]
    late_fluct = (float(np.std(ke_late) / np.mean(ke_late))
                  if len(ke_late) >= 2 and np.mean(ke_late) > 1e-15 else None)
    coherent_late_fluct = bool(late_fluct is not None and late_fluct < 0.05)

    l3_all_four_pass = bool(ke_nonzero and circ3d_nonzero and flux_nonzero and coherent_late_fluct)
    l3_partial = bool((ke_nonzero or circ3d_nonzero or flux_nonzero or coherent_late_fluct)
                       and not l3_all_four_pass)

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = True
    detected["localization"] = growing
    detected["spontaneous_motion"] = com_velocity_nonzero
    detected["circulation"] = circulation_nonzero

    # reached_level=3は全要件(KE≠0 AND 循環≠0 AND flux(Nu>1) AND coherent)を揃えた場合のみ。
    # 一部だけなら reached=2 のまま、l3_partial=True で明記する(前回round L4/L7の教訓)。
    reached = 2
    if l3_all_four_pass:
        reached = 3

    # coherentゲートは「構造が保たれたまま動くか」を問うもの——構造がまだ存在しない初期
    # ノイズ遷移(多モード競合で見かけ上coherenceが低い)ではなく、後期(飽和/単一ロールが
    # 支配的になった後)の値で判定する。ルートBの early/late 分割と同じ考え方。
    coherent_gate = bool(late_mean_coherence is not None and late_mean_coherence > 0.3)
    role = "E" if reached >= 3 else "F"
    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {
            "Ra": Ra, "Ra_c_theory": RA_C_THEORY, "Pr": p["Pr"], "grid": [NZ, NY, NX],
            "t_final": t_final, "dt": dt, "diverged": phys["diverged"],
            "w_rms_series": [round(x, 6) if np.isfinite(x) else None for x in w_rms_series],
            "circulation_proxy_series": [round(x, 6) for x in circ_series],
            "fitted_growth_rate": slope, "growing": growing,
            "circulation_nonzero": circulation_nonzero,
            "coherence_series": [round(c, 4) for c in coherence_series],
            "mean_coherence": round(mean_coherence, 4) if mean_coherence is not None else None,
            "early_phase_mean_coherence": round(early_mean_coherence, 4) if early_mean_coherence is not None else None,
            "late_phase_mean_coherence": round(late_mean_coherence, 4) if late_mean_coherence is not None else None,
            "coherent_gate_pass": coherent_gate,
            "ke_series": [round(x, 6) for x in ke_series],
            "circulation_3d_series": [round(x, 6) for x in circ3d_series],
            "nu_series": [round(x, 6) for x in nu_series],
            "nu_late_mean": round(nu_late_mean, 6),
            "late_fluct": round(late_fluct, 6) if late_fluct is not None else None,
            "l3_gate_full": {
                "com_velocity_nonzero": ke_nonzero,
                "circulation_nonzero": circ3d_nonzero,
                "flux_nonzero_Nu_gt_1": flux_nonzero,
                "coherent_late_fluct_lt_0p05": coherent_late_fluct,
                "all_four_pass": l3_all_four_pass,
            },
            "l3_partial": l3_partial,
            "l3_gate": {"com_velocity_nonzero": com_velocity_nonzero,
                       "circulation_nonzero": circulation_nonzero,
                       "coherent_while_growing": coherent_gate},
            "honest_limitation": ("axis順序バグ修正後は数値発散なく非線形飽和(w_rmsの伸びの"
                                  "頭打ち)まで到達する。coherentフラグは依頼書指定のlate_fluct"
                                  "(KEの後期相対std)<0.05で判定する——完全に時間的に一定な"
                                  "定常ロールに到達したかの直接指標。位相相関ベースのcoherence"
                                  "(early/late split)は補助的な第二の指標として併記する。"
                                  "Nu(=1+<w*T>)は後期(最終1/4区間)の時間平均で評価する。"),
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": role},
    }
    integrity = io.integrity_block(
        conservation_drift=0.0,  # Boussinesqは厳密な保存則なし(壁を通じた熱流入出があり得る、該当なし)
        resolutions_result={"%dx%dx%d" % (NZ, NY, NX): reached},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False,  # 静止+Tノイズのみ、対流セル/回転/循環は置かない
        gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "Ra=300 (<Ra_c=657.5, decisive frozen control)",
                        "result": "see companion room -- expect decay, no motion"}])
    checksum = io.checksum_of(rb.physical_slice(snapshots[-1]["w"], NZ))
    genesis_yaml = {
        "equations": "Boussinesq RB (free-slip, fixed-T), domain-doubling pseudo-spectral, "
                     "Ra_c(theory)=27*pi^4/4~657.51",
        "solver": "semi-implicit spectral (diffusion) + explicit (advection/buoyancy), "
                 "Chorin projection (Leray)",
        "dt": dt, "dx": None, "grid": [NZ, NY, NX], "boundary": "free-slip + fixed-T (z), periodic (x,y)",
        "params": p, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("L3-ルートA [%s] Ra=%.0f (Ra_c理論値=%.1f), t_final=%.2f, grid=%s。w_rms: %s。"
             "fitted_growth_rate=%s, growing=%s, circulation_nonzero=%s, mean_coherence=%s "
             "(early=%s, late=%s, coherent_gate=%s)。全要件ゲート: KE≠0=%s, 循環3D≠0=%s, "
             "Nu(後期)=%.4f(flux≠0=%s), late_fluct=%s(coherent=%s) -> all_four_pass=%s "
             "(l3_partial=%s)。diverged=%s。%s"
             % (notes_extra, Ra, RA_C_THEORY, t_final, (NZ, NY, NX),
                ["%.2e" % x for x in w_rms_series], slope, growing, circulation_nonzero,
                mean_coherence, early_mean_coherence, late_mean_coherence, coherent_gate,
                ke_nonzero, circ3d_nonzero, nu_late_mean, flux_nonzero, late_fluct,
                coherent_late_fluct, l3_all_four_pass, l3_partial,
                phys["diverged"], notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  Ra=%.0f KE_nonzero=%s circ3d_nonzero=%s Nu=%.4f late_fluct=%s "
          "all_four_pass=%s reached=%d"
          % (run_dir, Ra, ke_nonzero, circ3d_nonzero, nu_late_mean, late_fluct,
             l3_all_four_pass, reached))
    return snapshots, report


def main():
    print("=== L3-ルートA【最有力】Rayleigh-Bénard対流 full-3D ===")
    print("予言: Ra<657.5は静止(減衰)、Ra>657.5は自発的に成長->非線形飽和(coherent L3)。")
    print("全要件ゲート(KE≠0 AND 3D循環≠0 AND Nu>1 AND late_fluct<0.05)で判定。")

    print("\n[決定的対照: 亜臨界] Ra=300 (<Ra_c)")
    run_and_save("l3-route-a-subcritical-Ra300-seed0001", 300.0, seed=1, t_final=8.0,
                 notes_extra="決定的対照(亜臨界): 減衰、動かないことを期待(全要件不通過のはず)。"
                             "t_final=8まで延長し、真に定常ゼロに収束することを確認する。")

    print("\n[超臨界] Ra=1200 (>Ra_c)")
    run_and_save("l3-route-a-supercritical-Ra1200-seed0001", 1200.0, seed=1, t_final=8.0,
                 notes_extra="超臨界: coherentな成長->非線形飽和(定常ロール)を期待。"
                             "t_final=3では循環proxyがまだ緩やかに増加中だったため、t_final=8まで"
                             "延長し真の定常状態(late_fluct<0.05)への収束を確認する。")

    print("=== L3-ルートA done ===")


if __name__ == "__main__":
    main()
