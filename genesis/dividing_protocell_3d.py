#!/usr/bin/env python3
"""genesis/dividing_protocell_3d.py -- 依頼3: 分裂プロトセル（3D、保存二場 Zwicker 型）。

    d(u)/dt = M*lap(mu) + R,     mu = f'(u) - kappa*lap(u),   f'(u) = 2*W*u*(1-u)*(1-2u)
    d(v)/dt = D_v*lap(v) - R
    R = k_p*v*(1-chi(u)) - k_d*chi(u),   chi(u) = 0.5*(1+tanh((u-0.5)/w))

u（液滴材料、Cahn-Hilliard 相分離）と v（可溶材料、拡散のみ）は R によって相互変換され、
総量 (u+v) は解析的に厳密保存される（同じ R が両式に符号違いで現れるため）。数値法は G003 と
同じ半陰的スペクトル（CH の biharmonic 剛性を陰的に、反応項は陽的に）。

このファイルは2種類の run を提供する：
  - `run`（natural genesis）: u は希薄バイアス平均 + 微小ノイズ、v は一様。完成した液滴は
    置かない（依頼書 §A initial_state）。「自然発生」の主張はこちらのみ。
  - `run_droplet_probe`: 単一液滴（半径 R0）を背景に置いた決定的対照/診断実験。分裂閾値・
    分裂モードを機構として特定するための**プローブ**であり、「自然発生」の主張には使わない
    （docs/PHYSICS_INTEGRITY.md の役割区別と同じ考え方：足場付き診断 vs 純粋創発）。
"""

import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

MODEL_ID = "dividing_protocell_3d"

DEFAULTS = {"M": 1.0, "kappa": 1.0, "W": 1.0, "w": 0.15, "D_v": 2.0, "k_p": 0.03, "k_d": 0.04}


def _chi(u, w):
    return 0.5 * (1.0 + np.tanh((u - 0.5) / w))


def _reaction(u, v, p):
    c = _chi(u, p["w"])
    return p["k_p"] * v * (1.0 - c) - p["k_d"] * c


def make_initial(shape, rng, u_mean=0.15, u_noise=0.02, v_uniform=0.30):
    """一様希薄バイアス + 微小ノイズ（u）、一様（v）。液滴は入れない（NO seeded droplets）。"""
    u = u_mean + u_noise * rng.standard_normal(shape)
    v = np.full(shape, v_uniform, dtype=float)
    return u, v


def make_droplet_initial(shape, R0, dx=1.0, center=None, u_background=0.15, u_inside=1.0,
                          v_uniform=0.30, interface_width=2.0, n_droplets=1, drop_spacing=None,
                          rng=None):
    """決定的対照/診断用：半径 R0（物理単位）の単一液滴（またはテスト用に複数個、周期像との
    干渉を避けて離して配置）を tanh 界面プロファイルで置く。自然発生の主張には使わない
    （本モジュール docstring）。R0・drop_spacing は物理単位（dx を掛けたグリッド座標で評価）。
    """
    ndim = len(shape)
    grids = np.indices(shape).astype(float) * dx  # physical coordinates
    u = np.full(shape, u_background, dtype=float)
    physical_size = [s * dx for s in shape]
    if n_droplets == 1:
        centers = [center or [s / 2.0 for s in physical_size]]
    else:
        centers = []
        drop_spacing = drop_spacing or (min(physical_size) / (n_droplets + 1))
        for i in range(n_droplets):
            c = [s / 2.0 for s in physical_size]
            c[0] = drop_spacing * (i + 1)
            centers.append(c)
    for c in centers:
        r = np.sqrt(sum((grids[a] - c[a]) ** 2 for a in range(ndim)))
        profile = 0.5 * (1.0 - np.tanh((r - R0) / interface_width))
        u = np.maximum(u, u_background + (u_inside - u_background) * profile)
    v = np.full(shape, v_uniform, dtype=float)
    if rng is not None:
        u = u + 0.001 * rng.standard_normal(shape)
    return u, v


def step(u, v, dt, p, k2):
    R = _reaction(u, v, p)
    fprime = 2.0 * p["W"] * u * (1.0 - u) * (1.0 - 2.0 * u)
    fprime_hat = np.fft.fftn(fprime)
    R_hat = np.fft.fftn(R)
    uhat = np.fft.fftn(u)
    vhat = np.fft.fftn(v)

    denom_u = 1.0 + dt * p["M"] * p["kappa"] * k2 ** 2
    uhat_new = (uhat + dt * (-p["M"] * k2 * fprime_hat + R_hat)) / denom_u
    denom_v = 1.0 + dt * p["D_v"] * k2
    vhat_new = (vhat - dt * R_hat) / denom_v

    u_new = np.real(np.fft.ifftn(uhat_new))
    v_new = np.real(np.fft.ifftn(vhat_new))
    return u_new, v_new


def _run_core(u0, v0, shape, steps, dt, params, snapshot_every, dx=1.0):
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    _, k2 = k_grid(shape, spacing=dx)
    u, v = u0, v0
    mass0 = float(np.mean(u + v))
    snapshots = []
    diverged = False
    v_bg_series = []
    for t in range(steps):
        u, v = step(u, v, dt, p, k2)
        if not np.all(np.isfinite(u)) or not np.all(np.isfinite(v)):
            diverged = True
            break
        v_bg_series.append(float(np.mean(v)))
        if t % snapshot_every == 0 or t == steps - 1:
            snapshots.append({"step": t, "u": u.copy(), "v": v.copy()})
    mass1 = float(np.mean(u + v))
    phys = {"mass_drift": mass1 - mass0, "diverged": diverged, "v_background_series": v_bg_series,
            "v_background_final": v_bg_series[-1] if v_bg_series else float(np.mean(v))}
    if not snapshots:
        snapshots = [{"step": 0, "u": u.copy(), "v": v.copy()}]
    return snapshots, phys


def run(shape, steps, dt, seed, params=None, snapshot_every=None, u_mean=0.15, u_noise=0.02,
        v_uniform=0.30, dx=1.0):
    """natural genesis: t=0（希薄バイアス+ノイズ、完成した液滴なし）から中断なく発展させる。"""
    rng = np.random.default_rng(seed)
    u0, v0 = make_initial(shape, rng, u_mean=u_mean, u_noise=u_noise, v_uniform=v_uniform)
    snapshot_every = snapshot_every or max(1, steps // 20)
    return _run_core(u0, v0, shape, steps, dt, params, snapshot_every, dx=dx)


def run_droplet_probe(shape, R0, steps, dt, seed, params=None, snapshot_every=None,
                       u_background=0.15, v_uniform=0.30, n_droplets=1, dx=1.0):
    """決定的対照/診断: 半径 R0（物理単位）の単一（または複数の隔離された）液滴から発展させる。"""
    rng = np.random.default_rng(seed)
    u0, v0 = make_droplet_initial(shape, R0, dx=dx, u_background=u_background, v_uniform=v_uniform,
                                   n_droplets=n_droplets, rng=rng)
    snapshot_every = snapshot_every or max(1, steps // 30)
    return _run_core(u0, v0, shape, steps, dt, params, snapshot_every, dx=dx)


def droplet_components(u, thr=0.5):
    """u>thr の連結成分（液滴）を数え、各体積（voxel 数）を返す。"""
    count, labeled, sizes = diag.connected_components(u > thr)
    return count, sizes


# --- Room 依頼3 (分裂プロトセル) 実験計画 -----------------------------------------------

DT = 0.05
STEPS = 3000  # t_final = 150
N_HI = 64
N_LO_L64 = 48  # same physical domain L=64 as N_HI, coarser dx -- genuine resolution check
CONFIRMED_PARAMS = dict(DEFAULTS)  # 2D-confirmed: M=1,kappa=1,W=1,w=0.15,D_v=2,k_p=0.03,k_d=0.04


def _analyze_components(snapshots, dx=1.0, thr=0.5):
    times, counts, size_lists = [], [], []
    for s in snapshots:
        n, sizes = droplet_components(s["u"], thr=thr)
        times.append(s["step"] * DT)
        counts.append(n)
        size_lists.append(sorted((sizes * dx ** 3).tolist(), reverse=True))
    return times, counts, size_lists


def _classify_division(counts, initial_count, max_burst_ratio=2.0, max_burst_absolute=3):
    """count 系列から「分裂」候補を検出するが、多数の成分がほぼ同時に現れる核形成バーストは
    分裂と区別する（依頼書の背景セクションが警告する「核形成と分裂の交絡」）。バースト＝直前比
    max_burst_ratio 倍を超え、かつ絶対増分が max_burst_absolute を超える単一スナップショット間の
    増加。バーストが一度でも検出されたら raw な増加があっても divided=False とする。
    """
    final = counts[-1]
    maxc = max(counts)
    nucleation_burst_detected = False
    for i in range(1, len(counts)):
        prev, cur = counts[i - 1], counts[i]
        if cur <= prev:
            continue
        absolute = cur - prev
        ratio = float("inf") if prev == 0 else cur / prev
        if absolute > max_burst_absolute and ratio > max_burst_ratio:
            nucleation_burst_detected = True
    raw_increase = bool(maxc > initial_count)
    divided = bool(raw_increase and not nucleation_burst_detected)
    div_idx = next((i for i, c in enumerate(counts) if c > initial_count), None)
    return {"divided": divided, "division_snapshot_index": div_idx, "initial_count": initial_count,
            "final_count": final, "max_count": maxc,
            "nucleation_burst_detected": nucleation_burst_detected, "raw_count_increase": raw_increase}


def _far_from_center_components(u_field, shape, dx=1.0, thr=0.5, center=None, max_distance_factor=2.5,
                                 R0=None):
    """probe run 用：連結成分の重心が液滴の初期配置中心から離れすぎていないかを見る（分裂で
    できた娘液滴 vs 背景の別の場所で勝手に核形成した液滴、を空間的に区別する）。
    Returns (n_far, n_total) -- n_far>0 なら背景核形成の疑いがある成分が存在。
    """
    count, labeled, sizes = diag.connected_components(u_field > thr)
    if count == 0:
        return 0, 0
    physical_center = center or [s * dx / 2.0 for s in shape]
    coms = np.atleast_2d(ndimage.center_of_mass(np.ones_like(labeled), labeled,
                                                  index=np.arange(1, count + 1)))
    threshold = max_distance_factor * (R0 if R0 else 5.0)
    n_far = 0
    for com in coms:
        phys_com = [c * dx for c in com]
        dist = float(np.sqrt(sum((phys_com[a] - physical_center[a]) ** 2 for a in range(len(shape)))))
        if dist > threshold:
            n_far += 1
    return n_far, count


def _save_room(room_id, genesis_yaml, emergence_report, integrity, input_vs_output, notes):
    run_dir = io.write_results(room_id, genesis_yaml, emergence_report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  reached_level=%d role=%s" % (run_dir, emergence_report["reached_level"],
          emergence_report["purity"]["role"]))
    return run_dir


def run_natural_and_save(room_id, shape, seed, steps=STEPS, notes_extra=""):
    snapshots, phys = run(shape, steps, DT, seed, params=CONFIRMED_PARAMS)
    times, counts, size_lists = _analyze_components(snapshots)
    div = _classify_division(counts, counts[0])
    u_fields = [s["u"] for s in snapshots]
    _, growth_rate = diag.variance_growth(u_fields)
    _, prominence = diag.structure_factor_peak(u_fields[-1])

    v_bg = phys["v_background_series"]
    v_bg_bounded = bool(max(v_bg) < 1.0) if v_bg else True  # sanity ceiling; no runaway supersaturation

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = bool(growth_rate > 0 and prominence > 1.5)
    detected["localization"] = bool(div["max_count"] >= 1)
    detected["persistent_individuality"] = bool(div["final_count"] >= 1 and not phys["diverged"])
    # 核形成バーストは分裂ではない（依頼書「核形成と分裂の交絡」警告）: _classify_division が
    # nucleation_burst_detected の場合 divided=False にしているので、そのまま使う。
    detected["growth_division_inheritance"] = bool(div["divided"])
    reached = 0
    if detected["difference"]:
        reached = 1
    if detected["localization"]:
        reached = max(reached, 2)
    if detected["persistent_individuality"]:
        reached = max(reached, 4)
    if detected["growth_division_inheritance"]:
        reached = max(reached, 7)

    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"mass_drift": round(float(phys["mass_drift"]), 15),
                        "variance_growth": round(float(growth_rate), 6),
                        "structure_factor_prominence": round(float(prominence), 4),
                        "v_background_initial": round(v_bg[0], 4) if v_bg else None,
                        "v_background_final": round(phys["v_background_final"], 4),
                        "v_background_max": round(max(v_bg), 4) if v_bg else None,
                        "v_background_bounded": v_bg_bounded,
                        "droplet_count_final": div["final_count"],
                        "droplet_count_max": div["max_count"],
                        "spontaneous_nucleation_occurred": bool(div["max_count"] >= 1),
                        "nucleation_burst_detected": div["nucleation_burst_detected"],
                        "raw_count_increase_before_burst_filter": div["raw_count_increase"]},
        "purity": {"per_object_labels": False, "external_optimum": False,
                   "role": "E" if reached >= 1 else "F"},
    }
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"], resolutions_result={"%dx%dx%d" % shape: div["final_count"]},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False, gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "natural genesis (uniform+noise, NO seeded droplet)",
                        "result": "droplet_count_max=%d, nucleation_burst_detected=%s, "
                                  "v_background bounded=%s (max=%.3f)"
                                  % (div["max_count"], div["nucleation_burst_detected"], v_bg_bounded,
                                     max(v_bg) if v_bg else 0.0)}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]])
    genesis_yaml = {
        "equations": "du/dt=M*lap(mu)+R, mu=f'(u)-kappa*lap(u); dv/dt=D_v*lap(v)-R; "
                     "R=k_p*v*(1-chi(u))-k_d*chi(u)",
        "solver": "pseudo-spectral, semi-implicit CH (u) + implicit diffusion (v), explicit reaction",
        "dt": DT, "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": CONFIRMED_PARAMS, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼3 [検証] natural genesis（希薄バイアス+ノイズ、完成した液滴なし）, t_final=%.0f, "
             "grid=%s, seed=%d。droplet_count time series: %s。%s"
             % (steps * DT, shape, seed, list(zip([round(t) for t in times], counts)), notes_extra))
    _save_room(room_id, genesis_yaml, report, integrity, input_vs_output, notes)
    return snapshots, phys, report, times, counts


def run_probe_and_save(room_id, shape, R0, seed, steps=STEPS, dx=1.0, params=None, notes_extra=""):
    p = params or CONFIRMED_PARAMS
    snapshots, phys = run_droplet_probe(shape, R0, steps, DT, seed, params=p, dx=dx)
    times, counts, size_lists = _analyze_components(snapshots, dx=dx)
    # probe（単一液滴を置く）では「元の液滴付近から生じたか」という空間的基準の方が
    # 「一気に増えたか」という時間的基準（_classify_division の burst 検出、自然発生用）より
    # 直接的な交絡排除になる: raw_count_increase を使い、空間チェックのみで確度を判定する。
    div = _classify_division(counts, counts[0])
    n_far, n_total = _far_from_center_components(snapshots[-1]["u"], shape, dx=dx, R0=R0)
    background_nucleation_suspected = bool(n_far > 0)
    genuine_division = bool(div["raw_count_increase"] and not background_nucleation_suspected)
    mode = ("stable_no_division" if not genuine_division else
            "clean_two_way_division" if div["max_count"] == 2 else "multi_fragmentation")
    if div["raw_count_increase"] and background_nucleation_suspected:
        mode = "confounded_by_background_nucleation"

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["localization"] = True
    detected["persistent_individuality"] = bool(div["final_count"] >= 1 and not phys["diverged"])
    detected["growth_division_inheritance"] = genuine_division
    reached = 2
    if detected["persistent_individuality"]:
        reached = 4
    if detected["growth_division_inheritance"]:
        reached = 7

    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": False,  # プローブは液滴を初期条件に置く: t=0 一様ではない
        "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"mass_drift": round(float(phys["mass_drift"]), 15),
                        "R0": R0, "division_mode": mode, "division_snapshot_index": div["division_snapshot_index"],
                        "initial_count": div["initial_count"], "final_count": div["final_count"],
                        "max_count": div["max_count"], "nucleation_burst_detected": div["nucleation_burst_detected"],
                        "background_nucleation_suspected": background_nucleation_suspected,
                        "n_components_far_from_original_site": n_far,
                        "final_sizes_physical_volume": size_lists[-1][:8]},
        # 液滴は決定的対照/診断として置いている（自然発生の主張ではない）: role=S
        "purity": {"per_object_labels": True, "external_optimum": False, "role": "S"},
    }
    integrity = io.integrity_block(
        conservation_drift=phys["mass_drift"], resolutions_result={"%dx%dx%d(dx=%.3f)" % (shape[0], shape[1], shape[2], dx): div["final_count"]},
        seed_success={str(seed): reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # R0 はまさに調べたい量（分裂するか）そのもの: 明示的に申告
        gate_encodes_conclusion_causality=False, gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "R0 sweep (this room is one point in the sweep)",
                        "result": "division_mode=%s at R0=%.1f" % (mode, R0)}])
    checksum = io.checksum_of([snapshots[-1]["u"], snapshots[-1]["v"]])
    genesis_yaml = {
        "equations": "du/dt=M*lap(mu)+R, mu=f'(u)-kappa*lap(u); dv/dt=D_v*lap(v)-R; "
                     "R=k_p*v*(1-chi(u))-k_d*chi(u)",
        "solver": "pseudo-spectral, semi-implicit CH (u) + implicit diffusion (v), explicit reaction",
        "dt": DT, "dx": dx, "grid": list(shape), "boundary": "periodic",
        "params": p, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼3 [決定的対照/診断、自然発生の主張ではない] R0=%.1f probe, t_final=%.0f, grid=%s, "
             "dx=%.3f。droplet_count time series: %s。division_mode=%s。%s"
             % (R0, steps * DT, shape, dx, list(zip([round(t) for t in times], counts)), mode,
                notes_extra))
    _save_room(room_id, genesis_yaml, report, integrity, input_vs_output, notes)
    return snapshots, phys, report, times, counts, size_lists


def main():
    shape_hi = (N_HI, N_HI, N_HI)

    print("=== 依頼3 分裂プロトセル (3D, 保存二場 Zwicker) ===")
    print("[検証 1] natural genesis (背景核形成チェック), seeds=1,2,3, grid=%s" % (shape_hi,))
    for seed in (1, 2, 3):
        run_natural_and_save("dividing-protocell-3d-natural-seed%04d" % seed, shape_hi, seed=seed)

    print("[検証 2 / 核心 3] R0 sweep (決定的分裂テスト + 清潔ピンチ vs 多断片化), grid=%s" % (shape_hi,))
    R0_values = [5, 7, 9, 11, 13, 15, 17]
    for R0 in R0_values:
        room_id = "dividing-protocell-3d-probe-R0-%02d-seed0001" % R0
        run_probe_and_save(room_id, shape_hi, R0, seed=1,
                           notes_extra="R0 sweep の一点。臨界サイズ・分裂モード（清潔な1->2 vs "
                                       "多断片化）を特定するための決定的対照。")

    print("[解像度収束] R0=13 (臨界付近), 同一物理領域 L=%d, N=%d(dx=%.3f) vs N=%d(dx=1.0)"
          % (N_HI, N_LO_L64, N_HI / N_LO_L64, N_HI))
    run_probe_and_save("dividing-protocell-3d-probe-R0-13-res%d-seed0001" % N_LO_L64,
                       (N_LO_L64, N_LO_L64, N_LO_L64), 13, seed=1, dx=N_HI / N_LO_L64,
                       notes_extra="解像度収束チェック: 物理領域 L=%d を固定し N=%d vs N=%d(dx=1.0, "
                                   "上の R0=13 probe) で比較。" % (N_HI, N_LO_L64, N_HI))

    print("[核心 4] grow-divide サイクル: 大きめ超臨界液滴を長時間積分")
    run_probe_and_save("dividing-protocell-3d-probe-R0-17-longrun-seed0001", shape_hi, 17, seed=1,
                       steps=STEPS * 2,
                       notes_extra="[探索] grow-divide サイクル探索: t_final を2倍に延長し、分裂後の"
                                   "娘液滴が再成長・再分裂するかを追跡。")

    print("[探索 5] regime 地図: (k_p, k_d) を振る（R0=13, grid 48^3 で軽量に）")
    shape_regime = (48, 48, 48)
    regime_points = [
        ("kp_low", dict(CONFIRMED_PARAMS, k_p=0.015)),
        ("kp_high", dict(CONFIRMED_PARAMS, k_p=0.06)),
        ("kd_low", dict(CONFIRMED_PARAMS, k_d=0.02)),
        ("kd_high", dict(CONFIRMED_PARAMS, k_d=0.08)),
        ("Dv_low", dict(CONFIRMED_PARAMS, D_v=0.5)),
        ("Dv_high", dict(CONFIRMED_PARAMS, D_v=6.0)),
    ]
    for label, params in regime_points:
        room_id = "dividing-protocell-3d-regime-%s-R0-13-seed0001" % label
        run_probe_and_save(room_id, shape_regime, 13, seed=1, steps=2000, params=params,
                           notes_extra="[探索] regime 地図の一点: %s。2D確認点からのパラメータ変化。"
                                       % label)

    print("=== 依頼3 分裂プロトセル done ===")


if __name__ == "__main__":
    main()
