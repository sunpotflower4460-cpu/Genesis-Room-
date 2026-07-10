#!/usr/bin/env python3
"""genesis/_template.py — Genesis Room 雛形（2D 複素 Ginzburg-Landau クエンチの最小例）。

Codex はこの型に従って 3D ソルバ（Model H / CGL / 分裂プロトセル 等）を書く：
  1. genesis condition（fields, laws, initial=uniform+noise, params, seed）を dict で定義する。
  2. t=0 から時間発展させる（外部から完成形を置かない・途中で機能を足さない）。
  3. 途中で snapshots を記録する。
  4. common.emergence.compute_level_report で到達 Level を測定する。
  5. common.diagnostics / common.io.integrity_block で整合性（保存量ドリフト・複数 seed 再現性）を測定する。
  6. common.io.write_results で §B 形式（summary.json / manifest.json / figures / NOTES.md）に保存する。

物理は最小の 2D 例（複素 TDGL クエンチ、docs/EMERGENCE_LEVELS.md Level 1-2 の Kibble-Zurek 型創発）。
--- ここに Codex が 3D ソルバを入れる（このファイル自体は 3D を計算しない） ---
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import emergence, io  # noqa: E402

# 1. genesis condition -------------------------------------------------------
GENESIS = {
    "schema_version": 1,
    "model": "template_2d_cgl_quench",
    "equations": "d(psi)/dt = eps(t)*psi - |psi|^2 psi + D*lap(psi)  (complex TDGL / CGL quench)",
    "solver": "explicit Euler, periodic Laplacian (np.roll)",
    "dimension": 2,
    "grid": [96, 96],
    "dx": 1.0,
    "dt": 0.05,
    "boundary": "periodic",
    "params": {"D": 1.0, "eps_final": 1.0, "quench_start": 0.0, "quench_duration": 10.0},
    "initial_state": {"type": "uniform_plus_noise", "mean_amplitude": 0.0, "noise_amplitude": 1.0e-2},
    "seed": 0,
    "steps": 400,
    "snapshot_every": 25,
}


def _laplacian(z):
    """Periodic Laplacian over all axes (2D or 3D), unit spacing."""
    out = -2 * z.ndim * z
    for ax in range(z.ndim):
        out = out + np.roll(z, 1, ax) + np.roll(z, -1, ax)
    return out


def _eps_of_t(t, p):
    """Quench protocol: eps ramps from -eps_final to +eps_final (physical, defined ahead of time)."""
    q0, qd, ef = p["quench_start"], p["quench_duration"], p["eps_final"]
    frac = min(max((t - q0) / qd, 0.0), 1.0) if qd > 0 else 1.0
    return ef * (2.0 * frac - 1.0)


def _free_energy(psi, p):
    grad2 = sum(np.abs(np.roll(psi, -1, ax) - psi) ** 2 for ax in range(psi.ndim))
    return float(np.mean(0.25 * np.abs(psi) ** 4 + 0.5 * p["D"] * grad2))


def run(genesis=GENESIS):
    """genesis condition から t=0 で一様+微小ノイズを作り、中断なく発展させて snapshots を返す。
    完成した模様・欠陥・波長は一切初期条件に置かない（docs/EMERGENCE_LEVELS.md 「seeded 禁止」）。
    """
    rng = np.random.default_rng(genesis["seed"])
    shape = tuple(genesis["grid"])
    noise = genesis["initial_state"]["noise_amplitude"]
    psi = (noise * (rng.standard_normal(shape) + 1j * rng.standard_normal(shape))).astype(np.complex128)

    p = genesis["params"]
    dt = genesis["dt"]
    steps = genesis["steps"]
    snap_every = genesis["snapshot_every"]

    f0 = _free_energy(psi, p)
    snapshots = []
    for t in range(steps):
        eps = _eps_of_t(t * dt, p)
        psi = psi + dt * (eps * psi - (np.abs(psi) ** 2) * psi + p["D"] * _laplacian(psi))
        if t % snap_every == 0 or t == steps - 1:
            snapshots.append({"step": t, "field": psi.copy()})
    f1 = _free_energy(psi, p)
    return snapshots, {"conservation_drift": f1 - f0}


def run_with_seeds(seeds):
    """複数 seed で再現性を確認する対照（docs/PHYSICS_INTEGRITY.md §2「複数 seed で再現性を確認」）。"""
    results = {}
    for seed in seeds:
        g = dict(GENESIS, seed=seed)
        snapshots, phys = run(g)
        report = emergence.compute_level_report(snapshots, kind="cgl")
        results[seed] = {"report": report, "phys": phys, "final_field": snapshots[-1]["field"]}
    return results


def main():
    seeds = [0, 1, 2]
    per_seed = run_with_seeds(seeds)

    primary = per_seed[seeds[0]]
    report = primary["report"]

    resolutions_result = {"96x96": report["measured_by"]["structure_factor_peak_k"]}
    seed_success = {str(s): per_seed[s]["report"]["reached_level"] for s in seeds}
    nan_or_clip = any(not np.all(np.isfinite(per_seed[s]["final_field"])) for s in seeds)
    integrity = io.integrity_block(
        conservation_drift=primary["phys"]["conservation_drift"],
        resolutions_result=resolutions_result,
        seed_success=seed_success,
        nan_or_clip=nan_or_clip,
    )

    # 第8監査（docs/PHYSICS_INTEGRITY.md §6）: 求める量を初期条件やゲートに埋め込んでいないことの
    # 自己点検。この雛形は温度パラメータ eps を時間の関数として決めているだけで、欠陥数・波長など
    # 測定対象そのものはどこにも入力していない。
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False,
        gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False,
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "no_quench_control (eps held negative)",
                        "result": "documented expectation only in this template: "
                                  "no symmetry breaking, no defects -- Codex should actually run this "
                                  "control for real claims, per docs/PHYSICS_INTEGRITY.md §6."}],
    )

    checksum = io.checksum_of(primary["final_field"])
    genesis_manifest = dict(GENESIS, checksum=checksum, seeds=seeds, commit=None)

    run_id = "template-2d-cgl-seed%04d" % GENESIS["seed"]
    notes = ("2D 例（雛形）。Codex はこの型で 3D ソルバ（genesis/ 配下）を書き、"
             "同じ common.emergence / common.io で §B 出力を作ること。"
             "この雛形自体は 2D の候補発見に過ぎず、正式結果ではない（docs/DIMENSION_POLICY.md）。")
    run_dir = io.write_results(
        run_id=run_id,
        genesis_yaml=genesis_manifest,
        emergence_report=report,
        integrity=integrity,
        input_vs_output=input_vs_output,
        figures={},
        notes=notes,
    )
    print("reached_level=%s (role=%s) wrote %s"
          % (report["reached_level"], report["purity"]["role"], run_dir))
    return run_dir


if __name__ == "__main__":
    main()
