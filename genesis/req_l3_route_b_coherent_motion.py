#!/usr/bin/env python3
"""L3-ルートB: g001渦線のcoherent移動/再結合（乱流でないL3）。

--- 予言（実行前に登録） ---
g001（CGL、L2確定: ノイズ→対称性破れ→渦線、full-3D）の渦線ネットワークを、from-0（ノイズ+
クエンチ、速度/回転/循環は一切imprintしない——第8監査L3）で長時間発展させる。前回の知見
「CGLは動くがchurn（清潔でない）」を、単なる欠陥数の活性(_classify_defect_dynamicsの
active)でなく、真のcoherence（構造が保たれたまま動くか）で厳密に確認する。

coherenceの測定: 位相相関（phase correlation, FFT crosscorrelation of consecutive defect
mask snapshots）で、隣接スナップショット間の最良の剛体シフトとそのピーク値（0〜1、1に近い
ほど「パターンがそのまま並進した」ことを意味し、パターンが再編成/churnした場合はピークが
低く広がる）を計算する。これはobject-tracking的な「同一構造が動くか」の代理指標——個々の
渦芯にIDを振る完全なtrackingではないが、パターン全体の剛体並進との一致度を直接測る。

予言: 欠陥ネットワークが密な初期段階（クエンチ直後）は turbulent（低coherenceピーク、
churn）。coarsening後の疎な段階でも、CGLは分散があるため反応拡散的な"個体"のような
自己伝播（渦リングのcurvature駆動並進等）を持たず、依然として低coherenceのまま
（=「動くがchurn」を追認）と予想する——ただし探索的に両フェーズを実測し、正直に判定する。

決定的対照: (b,c)=(0,0)を長時間(t_final=400)発展させた「凍る」regime（既知:
defect_tail_mean=0.0、静止・循環0のL2止まり）vs (b,c)=(0.3,0.3)の「動く」regime
（既知: 短時間ではdefect_tail_mean=303.9・active=True、ただし長時間でのcoherenceは未測定）。

falsification: もしcoarsening後の疎な段階でcoherenceピークが高い（構造が保たれたまま並進
する）なら、CGLでも清潔なcoherent L3が部分的に成立する、と訂正して報告する。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import diagnostics as diag, emergence, io  # noqa: E402
from genesis import g001_cgl_3d as cgl  # noqa: E402

SHAPE = (48, 48, 48)
T_FINAL = 400.0  # 既知: (0,0)はこの長さで完全にextinct(0.0)に達する


def defect_mask(field, amplitude_threshold_frac=0.25):
    """3方向のプラケット巻きを合算した3D欠陥密度マップ(実数、位置情報を保持)。
    common.diagnostics.winding_defect_count と同じ検出ロジックだが、スカラー合計でなく
    空間分布(ND配列)を返す——coherence(位相相関)の計算に使う。"""
    theta = np.angle(field)
    amp = np.abs(field)
    thr = amplitude_threshold_frac * float(amp.max()) if amp.size else 0.0
    density = np.zeros(field.shape, dtype=float)
    for axis_a, axis_b in [(1, 2), (0, 1), (0, 2)]:
        w = diag._plaquette_winding(theta, axis_a, axis_b, amp, thr)
        density += np.abs(w)
    return density


def phase_correlation_shift(density0, density1):
    """隣接スナップショット間の最良の剛体シフトとcoherenceピーク値(0〜1)。
    ピークが高い=パターンがそのまま並進(coherent)。低い/広い=パターンが再編成(churn)。"""
    if density0.sum() < 1e-9 or density1.sum() < 1e-9:
        return None, 0.0
    F0 = np.fft.fftn(density0)
    F1 = np.fft.fftn(density1)
    R = F0 * np.conj(F1)
    denom = np.abs(R)
    denom[denom < 1e-12] = 1e-12
    R = R / denom
    corr = np.real(np.fft.ifftn(R))
    peak_idx = np.unravel_index(np.argmax(corr), corr.shape)
    peak_val = float(corr[peak_idx])
    shift = tuple(int(idx if idx <= s // 2 else idx - s) for idx, s in zip(peak_idx, corr.shape))
    return shift, peak_val


def run_and_track(room_id, params, seed=1, t_final=T_FINAL, notes_extra=""):
    p = dict(cgl.DEFAULTS)
    p.update(params)
    dt = cgl.stable_dt(3, p["b"])
    steps = int(round(t_final / dt))
    snapshot_every = max(1, steps // 60)  # coherence追跡のため十分細かく

    rng = np.random.default_rng(seed)
    A = cgl.make_initial(SHAPE, 0.01, rng)  # 一様近ゼロ+ノイズのみ、速度/回転/循環は入れない
    snapshots = []
    diverged = False
    for t in range(steps):
        A = cgl.step(A, dt, p)
        if not np.all(np.isfinite(A)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            snapshots.append({"step": t, "field": A.copy()})
    if not snapshots:
        snapshots = [{"step": 0, "field": A.copy()}]

    defect_series = [diag.winding_defect_count(s["field"]) for s in snapshots]
    densities = [defect_mask(s["field"]) for s in snapshots]
    times = [s["step"] * dt for s in snapshots]

    coherence_series = []
    for i in range(1, len(densities)):
        shift, peak = phase_correlation_shift(densities[i - 1], densities[i])
        coherence_series.append({"t": round(times[i], 1), "shift": shift, "peak_coherence": round(peak, 4),
                                 "n_defects": defect_series[i]})

    n_half = len(coherence_series) // 2
    early = coherence_series[:n_half] if n_half else []
    late = coherence_series[n_half:] if n_half else []
    early_mean_coh = float(np.mean([c["peak_coherence"] for c in early])) if early else None
    late_mean_coh = float(np.mean([c["peak_coherence"] for c in late])) if late else None

    dyn = cgl._classify_defect_dynamics(defect_series, dt=dt, snapshot_every=snapshot_every)

    # coherent L3の判定: 欠陥が活性(active, 消滅していない) かつ coherenceピークが高い(構造維持)
    COHERENCE_THRESHOLD = 0.5  # 剛体並進なら理論上1.0に近い、無相関なら1/N(ボクセル数)程度に低い
    late_coherent = bool(late_mean_coh is not None and late_mean_coh > COHERENCE_THRESHOLD)
    motion_active = bool(dyn["active"])
    com_velocity_nonzero = bool(any(c["shift"] != (0, 0, 0) for c in coherence_series if c["peak_coherence"] > 0.1))
    circulation_nonzero = bool(defect_series[-1] > 0)

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = True
    detected["localization"] = bool(defect_series[0] > 0)
    detected["circulation"] = circulation_nonzero
    detected["spontaneous_motion"] = com_velocity_nonzero and motion_active
    # coherent(清潔な)L3固有のフラグ: 構造が保たれたまま動く証拠(位相相関ピーク高)
    coherent_l3 = bool(detected["spontaneous_motion"] and late_coherent)

    reached = 2
    if detected["localization"]:
        reached = 2
    if detected["spontaneous_motion"] and detected["circulation"]:
        reached = 3  # ゲート通過(com_velocity!=0 AND circulation!=0)——ただしcoherentかは別途明記

    role = "E" if reached >= 3 else "F"
    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {
            "params": p, "t_final": t_final, "grid": list(SHAPE),
            "defect_dynamics": dyn,
            "early_phase_mean_coherence": round(early_mean_coh, 4) if early_mean_coh is not None else None,
            "late_phase_mean_coherence": round(late_mean_coh, 4) if late_mean_coh is not None else None,
            "coherence_threshold": COHERENCE_THRESHOLD,
            "late_phase_judged_coherent": late_coherent,
            "coherent_l3_gate": {"com_velocity_nonzero": com_velocity_nonzero,
                                 "circulation_nonzero": circulation_nonzero,
                                 "motion_active_not_extinct": motion_active,
                                 "structure_coherent_while_moving": late_coherent,
                                 "all_four_pass": coherent_l3},
            "coherence_series_sampled": coherence_series[::max(1, len(coherence_series) // 20)],
            "defect_count_final": defect_series[-1],
            "honest_note": ("coherenceは位相相関ピーク(隣接スナップショット間の最良剛体シフト"
                            "との一致度)で測る代理指標であり、個々の渦芯へのID付与による厳密な"
                            "object-trackingではない。ピークが低い/広い場合、パターンが剛体並進"
                            "でなく再編成(churn)したことを意味する。"),
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": role},
    }
    integrity = io.integrity_block(
        conservation_drift=0.0, resolutions_result={"48x48x48": defect_series[-1]},
        seed_success={str(seed): reached}, nan_or_clip=diverged)
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False,  # 一様近ゼロ+ノイズのみ、渦は一切imprintしない
        gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "(b,c)=(0,0) frozen control (decisive: known extinct at t_final=400)",
                        "result": "see companion room -- expect defect_tail_mean=0, no motion"}])
    checksum = io.checksum_of(snapshots[-1]["field"])
    genesis_yaml = {
        "equations": "dA/dt = A + (1+i*b)*lap(A) - (1+i*c)*|A|^2*A  (CGL, 法則固定)",
        "solver": "real-space finite-difference (periodic np.roll Laplacian), explicit Euler",
        "dt": dt, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
        "params": p, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("L3-ルートB [%s] (b,c)=(%.2f,%.2f), t_final=%.0f。defect: %d->%d (active=%s, "
             "still_coarsening=%s, extinct=%s)。coherence: early=%s late=%s "
             "(threshold=%.2f, late_coherent=%s)。coherent_l3_gate全通過=%s。%s"
             % (notes_extra, p["b"], p["c"], t_final, defect_series[0], defect_series[-1],
                dyn["active"], dyn["still_coarsening"], dyn["extinct"], early_mean_coh, late_mean_coh,
                COHERENCE_THRESHOLD, late_coherent, coherent_l3, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  defect %d->%d active=%s coherent(late)=%s reached=%d"
          % (run_dir, defect_series[0], defect_series[-1], dyn["active"], late_coherent, reached))
    return snapshots, report


def main():
    print("=== L3-ルートB: g001渦線のcoherent移動/再結合 ===")
    print("予言: 密な初期段階はturbulent(低coherence)。疎なcoarsening後も分散ありCGLは")
    print("      自己伝播機構を持たず低coherenceのまま(要実測)。")

    print("\n[決定的対照: 凍るregime] (b,c)=(0,0), t_final=400 (既知: 完全にextinct)")
    run_and_track("l3-route-b-frozen-control-b0c0-seed0001", {"b": 0.0, "c": 0.0}, seed=1,
                 notes_extra="決定的対照(凍る): defect活性ゼロ、循環0、動かないことを期待。")

    print("\n[動くregime候補] (b,c)=(0.3,0.3), t_final=400")
    run_and_track("l3-route-b-active-b0p3c0p3-seed0001", {"b": 0.3, "c": 0.3}, seed=1,
                 notes_extra="動くregime候補: 短時間ではactive=True(defect_tail_mean=303.9)。"
                             "長時間でのcoherence(構造保持)を実測する。")

    print("\n[前回round: BF不安定regime再検証] (b,c)=(1.052,-1.684), t_final=400")
    run_and_track("l3-route-b-bf-unstable-b1p052cm1p684-seed0001", {"b": 1.052, "c": -1.684}, seed=1,
                 notes_extra="前回round(法則クラスまたぎ探索)でsustained_multi_grain_fraction=0.88"
                             "(t_final=60, 見かけ上sustained-turbulent)だったregime。t_final=400まで"
                             "延長して真に持続するか再検証する。")

    print("=== L3-ルートB done ===")


if __name__ == "__main__":
    main()
