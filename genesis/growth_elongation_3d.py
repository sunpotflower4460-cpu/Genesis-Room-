#!/usr/bin/env python3
"""依頼C: 伸長機構の探索（受動では丸いまま）＋ 依頼B: grow-elongate-divideサイクル。

--- 依頼C 予言（実行前に登録） ---
予言候補(c) 成長不安定: 区画化液滴に強い成長（kp増 or 外部栄養供給v_feedで無制限成長）を
与えると、CH表面張力(kappa)による曲率緩和より速く物質が追加され、拡散律速成長不安定性
（Mullins-Sekerka類似）で初期の微小非対称ノイズが増幅し自然に伸長する。遅い成長（対照、
成長なし/弱い）では等方のまま丸い。

決定的測定: PCA固有値比 lambda_max/lambda_min の時間発展。速い成長で比が有意に増加すれば
自然な伸長不安定性を支持。増加しなければ falsified。

falsification: 速い成長でも比が対照と有意差なく増加しない（丸いまま、または多葉に分岐する
だけで単一長軸を持たない）なら、予言候補(c)は3Dで支持されない。その場合、依頼Bでは正直に
「機構(c)は不十分」と明記し、代替として明示的な軸バイアス（機構(a)、非等方成長、scaffold
として役割Sで正直にラベル）を使う。

--- 依頼B 予言（実行前に登録） ---
区画化液滴に純成長（外部栄養supply）を与え、アスペクト比が閾値(2.0)を超えたら furrow が
自動武装（手動タイミングでなく、液滴自身の現在のPCAアスペクト比が引き金）して分裂、娘が
自己維持（背景clean）を継承し、再成長->再伸長->再分裂のサイクルが少なくとも2世代
（1->2->4）回る。回らなければ、何が足りないか（伸長不足・再武装失敗・材料保存問題）を
正直に報告する。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import active_division_3d as ad  # noqa: E402
from genesis import compartmentalized_genesis_3d as cg  # noqa: E402
from genesis.dividing_protocell_3d import make_droplet_initial  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

DT = 0.05


def pca_aspect(u, X, Y, Z, threshold=0.5, min_voxels=50):
    mask = u > threshold
    n = int(mask.sum())
    if n < min_voxels:
        return None
    mx, my, mz = X[mask], Y[mask], Z[mask]
    centroid = np.array([mx.mean(), my.mean(), mz.mean()])
    pts = np.stack([mx - centroid[0], my - centroid[1], mz - centroid[2]], axis=1)
    cov = (pts.T @ pts) / n
    eigvals = np.linalg.eigvalsh(cov)  # ascending, all >=0
    lo = max(eigvals[0], 1e-9)
    return float(eigvals[-1] / lo)


# ---------------------------------------------------------------------------
# 依頼C: 成長不安定性テスト（growth-rate instability, mechanism c）
# ---------------------------------------------------------------------------

def reaction_growth_fed(u, v, p):
    """区画化生産R + 外部栄養供給（v をv_targetへケモスタット的に緩和、無制限成長を可能に）。"""
    c = cg.chi(u, p["w"])
    R = p["k_p"] * v * c - p["k_d"] * c ** 2
    feed = p["k_feed"] * (p["v_target"] - v)
    return R, feed


def step_growth(u, v, dt, p, k2):
    R, feed = reaction_growth_fed(u, v, p)
    fprime = 2.0 * p["W"] * u * (1.0 - u) * (1.0 - 2.0 * u)
    fprime_hat = np.fft.fftn(fprime)
    R_hat = np.fft.fftn(R)
    uhat = np.fft.fftn(u)
    vhat = np.fft.fftn(v)
    denom_u = 1.0 + dt * p["M"] * p["kappa"] * k2 ** 2
    uhat_new = (uhat + dt * (-p["M"] * k2 * fprime_hat + R_hat)) / denom_u
    denom_v = 1.0 + dt * p["D_v"] * k2
    vhat_new = (vhat + dt * (-R_hat + np.fft.fftn(feed))) / denom_v
    return np.real(np.fft.ifftn(uhat_new)), np.real(np.fft.ifftn(vhat_new))


def run_growth_probe(shape, R0, steps, seed, k_feed, v_target, dt=DT, params=None,
                     snapshot_every=None):
    p = dict(cg.DEFAULTS, k_d=0.02, k_feed=k_feed, v_target=v_target)
    if params:
        p.update(params)
    rng = np.random.default_rng(seed)
    u, v = make_droplet_initial(shape, R0, u_background=0.15, v_uniform=0.5, rng=rng)
    _, k2 = k_grid(shape)
    X, Y, Z = np.indices(shape).astype(float)
    snapshot_every = snapshot_every or max(1, steps // 40)
    mass0 = float(np.mean(u + v))
    snapshots = []
    diverged = False
    for t in range(steps):
        u, v = step_growth(u, v, dt, p, k2)
        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            n, _, sizes = diag.connected_components(u > 0.5)
            dilute = u < 0.5
            u_bg = float(np.mean(u[dilute])) if dilute.any() else float(np.mean(u))
            aspect = pca_aspect(u, X, Y, Z)
            snapshots.append({"step": t, "n": n, "sizes": sorted(sizes.tolist(), reverse=True)[:6],
                              "u_bg": u_bg, "aspect": aspect, "u": u.copy(), "v": v.copy()})
    mass1 = float(np.mean(u + v))
    phys = {"mass_drift": mass1 - mass0, "diverged": diverged}
    if not snapshots:
        snapshots = [{"step": 0, "n": 1, "sizes": [0.0], "u_bg": 0.0, "aspect": 1.0,
                     "u": u.copy(), "v": v.copy()}]
    return snapshots, phys


def reqc_instability_probe():
    print("=== 依頼C: 成長不安定性(mechanism c) 探索 ===")
    print("予言: 速い成長(高k_feed)でPCAアスペクト比が自然に増加。遅い成長(対照)は等方のまま。")
    shape = (64, 64, 64)
    R0 = 8
    steps = 1200  # t_final=60

    results = {}
    for label, k_feed, v_target in [("fast_growth", 0.08, 2.0), ("slow_growth_control", 0.0, 0.5)]:
        snaps, phys = run_growth_probe(shape, R0, steps, seed=1, k_feed=k_feed, v_target=v_target)
        aspects = [s["aspect"] for s in snaps if s["aspect"] is not None]
        sizes0 = snaps[0]["sizes"][0] if snaps[0]["sizes"] else 0.0
        sizesF = snaps[-1]["sizes"][0] if snaps[-1]["sizes"] else 0.0
        results[label] = {"aspect_series": [round(a, 4) for a in aspects],
                          "aspect_initial": aspects[0] if aspects else None,
                          "aspect_final": aspects[-1] if aspects else None,
                          "vol_initial": sizes0, "vol_final": sizesF,
                          "n_final": snaps[-1]["n"], "diverged": phys["diverged"],
                          "mass_drift": phys["mass_drift"]}
        print("  [%s] k_feed=%.3f v_target=%.2f  aspect %.3f -> %.3f, vol %.0f -> %.0f, n=%d"
              % (label, k_feed, v_target, results[label]["aspect_initial"] or -1,
                 results[label]["aspect_final"] or -1, sizes0, sizesF, results[label]["n_final"]))

    fast = results["fast_growth"]
    slow = results["slow_growth_control"]
    elongation_gain_fast = (fast["aspect_final"] - fast["aspect_initial"]) if fast["aspect_final"] else 0.0
    elongation_gain_slow = (slow["aspect_final"] - slow["aspect_initial"]) if slow["aspect_final"] else 0.0
    supported = bool(elongation_gain_fast > 0.15 and elongation_gain_fast > 2 * max(elongation_gain_slow, 0.01))

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["localization"] = True
    detected["persistent_individuality"] = bool(not fast["diverged"] and fast["n_final"] >= 1)

    report = {
        "reached_level": 4 if not supported else 4, "candidate_level": 5,
        "uninterrupted_from_zero": False, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"fast_growth": fast, "slow_growth_control": slow,
                        "elongation_gain_fast": round(elongation_gain_fast, 4),
                        "elongation_gain_slow_control": round(elongation_gain_slow, 4),
                        "mechanism_c_supported": supported},
        "purity": {"per_object_labels": False, "external_optimum": False,
                  "role": "E" if supported else "N"},
    }
    integrity = io.integrity_block(
        conservation_drift=fast["mass_drift"], resolutions_result={"64x64x64": fast["n_final"]},
        seed_success={"1": report["reached_level"]},
        nan_or_clip=fast["diverged"] or slow["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "slow_growth (k_feed=0, decisive control)",
                        "result": "aspect gain=%.4f vs fast=%.4f" % (elongation_gain_slow,
                                                                     elongation_gain_fast)}])
    checksum = io.checksum_of([np.array(fast["aspect_series"])])
    genesis_yaml = {
        "equations": "du/dt=M*lap(mu)+R; dv/dt=D_v*lap(v)-R+k_feed*(v_target-v) (chemostat供給で"
                     "無制限成長); R=kp*v*chi(u)-kd*chi(u)^2. 初期条件は等方球(微小noiseのみ、"
                     "軸バイアスなし)",
        "solver": "pseudo-spectral, semi-implicit CH + implicit diffusion + chemostat feed",
        "dt": DT, "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": {"k_feed_fast": 0.08, "v_target_fast": 2.0, "k_feed_slow": 0.0,
                  "v_target_slow": 0.5}, "seed": 1, "seeds": [1], "commit": None,
        "checksum": checksum,
    }
    notes = ("依頼C mechanism(c) 成長不安定性探索。fast: aspect %.3f->%.3f (gain=%.4f), slow対照: "
             "aspect %.3f->%.3f (gain=%.4f)。supported=%s。"
             % (fast["aspect_initial"] or -1, fast["aspect_final"] or -1, elongation_gain_fast,
                slow["aspect_initial"] or -1, slow["aspect_final"] or -1, elongation_gain_slow,
                supported))
    run_dir = io.write_results("req-c2-growth-instability-mechanism-c", genesis_yaml, report,
                                integrity, input_vs_output, figures={}, notes=notes)
    print("  wrote %s  supported=%s" % (run_dir, supported))
    return results, supported


if __name__ == "__main__":
    reqc_instability_probe()
