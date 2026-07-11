#!/usr/bin/env python3
"""genesis/g001_cgl_3d.py -- 依頼2: G001 複素 Ginzburg-Landau 欠陥形成 (3D)。

    dA/dt = A + (1 + i*b) lap(A) - (1 + i*c) |A|^2 A          (complex Ginzburg-Landau)

完成した渦線/ループは一切初期条件に置かない：A は一様近ゼロ + 微小ノイズから開始する
(requests/...md 依頼2 §A initial_state)。(b,c)=(0,0) は既に supercritical（+A の線形成長項）
なのでクエンチ・スケジュールは不要 -- t=0 から一定係数で積分する。

数値法：周期実空間 Laplacian（np.roll、次元非依存）、explicit Euler。(1+ib) 項の陽的離散化は
b が大きいほど硬くなる（拡散固有値の虚部が大きくなるため）ので、`stable_dt` で b から安全な dt を
計算する。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import emergence, io, diagnostics as diag  # noqa: E402

MODEL_ID = "g001_cgl_3d"

DEFAULTS = {"b": 0.0, "c": 0.0}


def stable_dt(ndim, b, dx=1.0, safety=0.4, dt_min=0.002, dt_max=0.05):
    """Explicit-Euler stability bound for dA/dt ⊃ (1+ib)*lap(A): at the most negative periodic-stencil
    eigenvalue (-4*ndim/dx^2), the complex growth rate is g = 1 - 4*ndim/dx^2 + i*b*(-4*ndim/dx^2);
    |1+dt*g| <= 1 requires dt <= -2*Re(g)/|g|^2. Larger |b| tightens this a lot (imaginary part grows).
    """
    lap_max = 4.0 * ndim / dx ** 2
    re_g = 1.0 - lap_max
    im_g = b * lap_max
    denom = re_g ** 2 + im_g ** 2
    dt = safety * (-2.0 * re_g) / denom if denom > 0 else dt_max
    return float(min(max(dt, dt_min), dt_max))


def make_initial(shape, noise_amplitude, rng):
    """一様近ゼロ + 微小ノイズ。渦線/ループは入れない（NO seeded vortex lines）。"""
    return (noise_amplitude * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))).astype(np.complex128)


def _laplacian(z):
    out = -2 * z.ndim * z
    for ax in range(z.ndim):
        out = out + np.roll(z, 1, ax) + np.roll(z, -1, ax)
    return out


def step(A, dt, p):
    lap = _laplacian(A)
    return A + dt * (A + (1 + 1j * p["b"]) * lap - (1 + 1j * p["c"]) * (np.abs(A) ** 2) * A)


def run(shape, t_final, seed, params=None, snapshot_every_frac=1.0 / 20, noise_amplitude=0.01, dt=None):
    """t=0（一様近ゼロ+微小ノイズ）から中断なく発展させ、snapshots を返す。dt は None なら b から
    自動選択（stable_dt）。"""
    p = dict(DEFAULTS)
    if params:
        p.update(params)
    ndim = len(shape)
    dt = dt if dt is not None else stable_dt(ndim, p["b"])
    steps = int(round(t_final / dt))
    snapshot_every = max(1, int(steps * snapshot_every_frac))

    rng = np.random.default_rng(seed)
    A = make_initial(shape, noise_amplitude, rng)

    snapshots = []
    diverged = False
    for t in range(steps):
        A = step(A, dt, p)
        if not np.all(np.isfinite(A)):
            diverged = True
            break
        if t % snapshot_every == 0 or t == steps - 1:
            snapshots.append({"step": t, "field": A.copy()})
    if not snapshots:
        snapshots = [{"step": 0, "field": A.copy()}]
    phys = {"diverged": diverged, "dt_used": dt, "steps": steps}
    return snapshots, phys


# --- Room 依頼2 (G001) 実験計画: 検証 + (b,c) 相図探索 + 3D トポロジー + 決定的対照 --------

T_FINAL = 100.0
N_HI = 48
N_LO = 32
NOISE_AMPLITUDE = 0.01


def _classify_defect_dynamics(defect_series, dt, snapshot_every):
    """欠陥数の時系列から動態を分類する。単純な tail の変動係数(cv)は、まだ粗大化/対消滅が
    進行中で単調減少している区間でも高い値を示してしまう（減少トレンドを"変動"と誤認する）ため
    使わない。代わりに tail を線形デトレンドし、(a) トレンドが tail 平均に対してまだ大きい
    （＝定常に達していない＝coarsening_still_in_progress）か、(b) トレンド除去後の残差変動が
    有意（＝定常状態での持続的な生成・消滅＝active）かを区別する。

    Returns dict with keys: tail_mean, still_coarsening, active, extinct.
    """
    tail = np.array(defect_series[len(defect_series) // 2:], dtype=float)
    tail_mean = float(np.mean(tail)) if len(tail) else 0.0
    if tail_mean < 1.0:
        return {"tail_mean": tail_mean, "still_coarsening": False, "active": False, "extinct": True,
                "trend_relative": 0.0, "residual_relative_std": 0.0}
    t_idx = np.arange(len(tail))
    slope, intercept = np.polyfit(t_idx, tail, 1)
    trend_relative = abs(slope * len(tail)) / tail_mean  # net fractional change in the tail window
    residual = tail - (slope * t_idx + intercept)
    residual_relative_std = float(np.std(residual)) / tail_mean
    still_coarsening = bool(trend_relative > 0.15)
    active = bool((not still_coarsening) and residual_relative_std > 0.05)
    return {"tail_mean": tail_mean, "still_coarsening": still_coarsening, "active": active,
            "extinct": False, "trend_relative": round(float(trend_relative), 4),
            "residual_relative_std": round(residual_relative_std, 4)}


def _level_report_cgl(snapshots):
    """compute_level_report の Lv1/Lv2 はそのまま使える（複素場・位相巻き）。CGL には陽な流速場が
    ないため、Lv3（自発運動/相互作用）は欠陥数の時系列の動態で代理測定する：`_classify_defect_dynamics`
    参照。まだ粗大化中（still_coarsening）の場合は「凍結」とも「活性」とも断定しない
    （measured_by に inconclusive として記録、reached_level は Lv3 に進めない）。
    """
    report = emergence.compute_level_report(snapshots, kind="cgl")
    defect_series = [diag.winding_defect_count(s["field"]) for s in snapshots]
    dyn = _classify_defect_dynamics(defect_series, dt=None, snapshot_every=None)
    active = bool(report["detected"]["localization"] and dyn["active"])

    report["detected"]["spontaneous_motion"] = active
    report["detected"]["circulation"] = active  # CGL に速度場はない: 同じ欠陥動態を証拠として共有
    report["measured_by"]["defect_count_series"] = defect_series
    report["measured_by"]["defect_tail_mean"] = round(dyn["tail_mean"], 3)
    report["measured_by"]["defect_dynamics_still_coarsening"] = dyn["still_coarsening"]
    report["measured_by"]["defect_dynamics_extinct"] = dyn["extinct"]
    report["measured_by"]["defect_dynamics_trend_relative"] = dyn["trend_relative"]
    report["measured_by"]["defect_dynamics_residual_relative_std"] = dyn["residual_relative_std"]

    reached = report["reached_level"]
    if active:
        reached = max(reached, 3)
    report["reached_level"] = reached
    report["candidate_level"] = min(reached + 1, 8)
    role = "E" if reached >= 1 else "F"
    report["purity"] = {"per_object_labels": False, "external_optimum": False, "role": role}
    return report


def _run_and_save(room_id, shape, params, seed, notes_extra="", noise_amplitude=NOISE_AMPLITUDE,
                   t_final=T_FINAL):
    snapshots, phys = run(shape, t_final, seed, params=params, noise_amplitude=noise_amplitude)
    report = _level_report_cgl(snapshots)

    resolutions_result = {"%dx%dx%d" % shape: report["measured_by"]["defect_tail_mean"]}
    integrity = io.integrity_block(
        conservation_drift=0.0,  # CGL は保存則を持たない（散逸系）; drift は該当なしとして 0 を記録
        resolutions_result=resolutions_result,
        seed_success={str(seed): report["reached_level"]},
        nan_or_clip=phys["diverged"],
    )
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False,
        gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "(b,c)=(0,0) baseline (relaxational/frozen limit)",
                        "result": "reference point for comparing activity across the (b,c) sweep"}],
    )
    checksum = io.checksum_of(snapshots[-1]["field"])
    genesis_yaml = {
        "equations": "dA/dt = A + (1+i*b)*lap(A) - (1+i*c)*|A|^2*A  (complex Ginzburg-Landau)",
        "solver": "real-space finite-difference (periodic np.roll Laplacian), explicit Euler, "
                  "dt auto-selected from b via stable_dt()",
        "dt": phys["dt_used"], "dx": 1.0, "grid": list(shape), "boundary": "periodic",
        "params": params, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("依頼2 [検証/探索]: (b,c)=%s, t_final=%.0f (steps=%d, dt=%.4f), grid=%s。"
             "3D の渦線再結合・結び目/絡み合いは common/diagnostics.winding_defect_count の"
             "貫通面カウント（線本数の下限プロキシ、docstring 記載）のみで評価しており、"
             "明示的な線追跡・結び目検出は未実装（未探索、依頼書の「足したら明記」に対応）。%s"
             % (params, t_final, phys["steps"], phys["dt_used"], shape, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output, figures={},
                                notes=notes)
    print("  wrote %s  reached_level=%d role=%s active=%s coarsening=%s extinct=%s "
          "defect_tail_mean=%.1f trend=%.3f resid=%.3f diverged=%s"
          % (run_dir, report["reached_level"], report["purity"]["role"],
             report["detected"]["spontaneous_motion"], report["measured_by"]["defect_dynamics_still_coarsening"],
             report["measured_by"]["defect_dynamics_extinct"], report["measured_by"]["defect_tail_mean"],
             report["measured_by"]["defect_dynamics_trend_relative"],
             report["measured_by"]["defect_dynamics_residual_relative_std"], phys["diverged"]))
    return snapshots, phys, report


def main():
    shape_hi = (N_HI, N_HI, N_HI)
    shape_lo = (N_LO, N_LO, N_LO)

    print("=== 依頼2 G001 CGL 欠陥形成 (3D) ===")
    print("[primary] (b,c)=(0,0), seeds=1,2,3, grid=%s" % (shape_hi,))
    for seed in (1, 2, 3):
        _run_and_save("g001-cgl-3d-b0c0-seed%04d" % seed, shape_hi, {"b": 0.0, "c": 0.0}, seed=seed)

    print("[resolution] (b,c)=(0,0), grid=%s vs %s, seed=1" % (shape_lo, shape_hi))
    _run_and_save("g001-cgl-3d-b0c0-res32-seed0001", shape_lo, {"b": 0.0, "c": 0.0}, seed=1,
                  notes_extra="解像度収束チェック（32^3 vs 48^3、同一 spacing=1.0 -- 注意: N のみ"
                              "変えているため物理箱サイズも変わっており、厳密な意味での「同一領域を"
                              "細かくする」解像度収束ではなく有限サイズ依存チェックに近い。")

    print("[long] (b,c)=(0,0), t_final=400 で本当に凍結するか（tail がまだ粗大化中でないか確認）")
    _run_and_save("g001-cgl-3d-b0c0-long-seed0001", shape_hi, {"b": 0.0, "c": 0.0}, seed=1,
                  t_final=400.0,
                  notes_extra="[検証] 標準 t_final=100 の run では tail がまだ単調減少中で「凍結」"
                              "と断定できなかったため、t_final=400 まで延長して本当に定常/凍結に"
                              "達するか確認する。")

    print("[explore] (b,c) 相図")
    bc_points = [(0.3, 0.3), (0.5, 0.5), (1.0, 0.5), (1.0, -1.5), (1.5, -1.0)]
    for b, c in bc_points:
        bf = 1 + b * c
        room_id = "g001-cgl-3d-b%s-c%s-seed0001" % (str(b).replace(".", "p").replace("-", "m"),
                                                       str(c).replace(".", "p").replace("-", "m"))
        _run_and_save(room_id, shape_hi, {"b": b, "c": c}, seed=1,
                      notes_extra="[探索] (b,c)=(%.1f,%.1f), 1+bc=%.2f (%s)。" %
                                  (b, c, bf, "Benjamin-Feir 不安定域" if bf < 0 else "安定域"))

    print("[explore] noise_amplitude 感度 (b,c)=(0,0)")
    for amp in (0.001, 0.1):
        room_id = "g001-cgl-3d-b0c0-noise%s-seed0001" % str(amp).replace(".", "p")
        _run_and_save(room_id, shape_hi, {"b": 0.0, "c": 0.0}, seed=1, noise_amplitude=amp,
                      notes_extra="[探索] noise_amplitude=%.3g（既定 0.01 との比較、KZ 的な"
                                  "初期条件依存の代理チェック。真のクエンチ速度スキャンではない点に注意）。"
                                  % amp)

    print("=== 依頼2 G001 done ===")


if __name__ == "__main__":
    main()
