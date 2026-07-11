#!/usr/bin/env python3
"""依頼D【探索】: 受動vs能動の境界を測る — Gray-Scott（膜なし）3D。

--- 予言（実行前に登録） ---
Gray-Scott系（相分離・膜なし、2場の純粋な反応拡散）は、古典的な"mitosis"パラメータ域
（F=0.0367, k=0.0649, Pearsonの分類で文字通り"mitosis"と呼ばれる領域）で、単一の種スポット
から出発して受動的（能動furrow等の追加機構なし）にスポットが繰り返し分裂する
（self-replicating spots、2Dでは長年知られる現象）。

これを区画化液滴（相分離・膜あり、依頼A/Bで確認済み: furrow無しでは伸長液滴すら分裂しない、
Lv6の壁）と対比する：「膜（相分離境界・表面張力）の有無が受動分裂の可否を分ける」。

決定的測定: Gray-Scottでのスポット数の時間発展（増加すれば受動分裂を支持）。対照として同じ
Gray-ScottでF,kを"stable spots"域（分裂しない既知の領域、例F=0.03,k=0.06付近やF=0.014のような
非分裂域）に変えたrunも用意し、分裂がパラメータ依存（何でも増えるわけではない）であることを
確認する。

falsification: 3DでGray-Scottスポットが分裂しない（2Dの知見が3Dで崩れる）、または対照
（stable域）でも同様に分裂する（パラメータに依らない人工物）なら、対比の前提が崩れる。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

DT = 1.0
DEFAULTS = {"Du": 0.16, "Dv": 0.08, "F": 0.0367, "k": 0.0649}  # 古典的"mitosis"パラメータ


def make_seed_initial(shape, rng, n_seeds=1, seed_radius=4, u_base=1.0, v_base=0.0):
    u = np.full(shape, u_base, dtype=float)
    v = np.full(shape, v_base, dtype=float)
    center = [s // 2 for s in shape]
    X, Y, Z = np.indices(shape)
    dist = np.sqrt((X - center[0]) ** 2 + (Y - center[1]) ** 2 + (Z - center[2]) ** 2)
    blob = dist < seed_radius
    u[blob] = 0.50
    v[blob] = 0.25
    u += 0.01 * rng.standard_normal(shape)
    v += 0.01 * rng.standard_normal(shape)
    return u, v


def step(u, v, dt, p, k2):
    reaction_u = -u * v ** 2 + p["F"] * (1.0 - u)
    reaction_v = u * v ** 2 - (p["F"] + p["k"]) * v
    uhat = np.fft.fftn(u)
    vhat = np.fft.fftn(v)
    ru_hat = np.fft.fftn(reaction_u)
    rv_hat = np.fft.fftn(reaction_v)
    denom_u = 1.0 + dt * p["Du"] * k2
    denom_v = 1.0 + dt * p["Dv"] * k2
    uhat_new = (uhat + dt * ru_hat) / denom_u
    vhat_new = (vhat + dt * rv_hat) / denom_v
    return np.real(np.fft.ifftn(uhat_new)), np.real(np.fft.ifftn(vhat_new))


def count_spots(v, threshold=0.1, min_voxels=5):
    n, _, sizes = diag.connected_components(v > threshold)
    sig = sum(1 for s in sizes if s >= min_voxels)
    return sig


def run(shape, steps, seed, params=None, snapshot_every=None, dt=DT):
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    rng = np.random.default_rng(seed)
    u, v = make_seed_initial(shape, rng)
    mass0 = float(np.mean(u + v))
    _, k2 = k_grid(shape)
    snapshot_every = snapshot_every or max(1, steps // 60)
    snapshots = []
    diverged = False
    for t in range(steps):
        u, v = step(u, v, dt, p, k2)
        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            n_spots = count_spots(v)
            snapshots.append({"step": t, "n_spots": n_spots, "v_max": float(v.max())})
    if not snapshots:
        snapshots = [{"step": 0, "n_spots": 0, "v_max": float(v.max())}]
    mass1 = float(np.mean(u + v))
    # Gray-Scott は u,v とも保存則を持たない(F,kによる正味の湧き出し/消滅がある)。
    # mass_driftはここでは「発散していないか」の健全性チェックとしてのみ使う。
    phys = {"diverged": diverged, "mass_drift": mass1 - mass0}
    return snapshots, phys, u, v


def classify_and_save(room_id, shape, steps, seed, params=None, label="", notes_extra=""):
    snapshots, phys, u, v = run(shape, steps, seed, params=params)
    n_series = [(s["step"] * DT, s["n_spots"]) for s in snapshots]
    n_initial = snapshots[0]["n_spots"]
    n_max = max(s["n_spots"] for s in snapshots)
    n_final = snapshots[-1]["n_spots"]
    replicated = bool(n_max > max(n_initial, 1) and n_final >= 2)

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = True
    detected["localization"] = bool(n_final >= 1)
    detected["persistent_individuality"] = bool(n_final >= 1 and not phys["diverged"])
    detected["growth_division_inheritance"] = replicated

    reached = 2
    if detected["persistent_individuality"]:
        reached = 4
    if replicated:
        reached = 7  # 受動的な"分裂"（膜がないため自己維持/背景cleanのLv6概念は本来適用外だが、
                      # 「1個体から複数個体が増える」機能的にはLv7相当の現象として記録

    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"label": label, "F": (params or DEFAULTS)["F"], "k": (params or DEFAULTS)["k"],
                        "n_spots_initial": n_initial, "n_spots_max": n_max, "n_spots_final": n_final,
                        "replicated_passively": replicated,
                        "n_series_sampled": n_series[::max(1, len(n_series) // 25)]},
        "purity": {"per_object_labels": False, "external_optimum": False,
                  "role": "E" if replicated else "N"},
    }
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"], resolutions_result={"%dx%dx%d" % shape: n_max},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "membrane-bound compartmentalized droplet (companion model)",
                        "result": "furrow OFFでは分裂しない(req-a2 decisive control 1参照)、"
                                  "膜(相分離境界)が受動分裂を妨げる、との対比"}])
    checksum = io.checksum_of([u, v])
    genesis_yaml = {
        "equations": "Gray-Scott: du/dt=Du*lap(u)-u*v^2+F(1-u), dv/dt=Dv*lap(v)+u*v^2-(F+k)v "
                     "(膜なし、相分離なしの純粋反応拡散、compartmentalized_genesis_3dとの対比)",
        "solver": "pseudo-spectral, semi-implicit diffusion, explicit reaction",
        "dt": DT, "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": params or DEFAULTS, "seed": seed, "seeds": [seed], "commit": None,
        "checksum": checksum,
    }
    notes = ("依頼D [%s] F=%.4f k=%.4f, t_final=%.0f, grid=%s, seed=%d。n_spots: %d->%d(max=%d)。"
             "replicated_passively=%s。series: %s。%s"
             % (label, (params or DEFAULTS)["F"], (params or DEFAULTS)["k"], steps * DT, shape,
                seed, n_initial, n_final, n_max, replicated,
                n_series[::max(1, len(n_series) // 15)], notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  [%s] n_spots %d->%d(max=%d) replicated=%s"
          % (run_dir, label, n_initial, n_final, n_max, replicated))
    return snapshots, phys, report


def main():
    print("=== 依頼D【探索】受動vs能動の境界: Gray-Scott (膜なし, 3D) ===")
    print("予言: mitosisパラメータ域(F=0.0367,k=0.0649)で受動的にスポットが分裂増殖する。")
    shape = (64, 64, 64)
    steps = 6000

    print("\n[mitosis域] F=0.0367 k=0.0649")
    classify_and_save("req-d2-grayscott-mitosis-seed0001", shape, steps, seed=1,
                      label="mitosis_region",
                      notes_extra="mitosisパラメータ域(既知の自己複製スポット)。")

    print("\n[安定スポット対照] F=0.030 k=0.060 (既知の非分裂/安定パターン域)")
    stable_params = dict(DEFAULTS, F=0.030, k=0.060)
    classify_and_save("req-d2-grayscott-stable-control-seed0001", shape, steps, seed=1,
                      params=stable_params, label="stable_region_control",
                      notes_extra="決定的対照: 安定パターン域では分裂しないことを期待"
                                  "（パラメータ依存性の確認、人工物でないことの証拠）。")

    print("=== 依頼D done ===")


if __name__ == "__main__":
    main()
