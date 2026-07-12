#!/usr/bin/env python3
"""確認B: 完全なL7を3D tagged Gray-Scottで（分裂＋遺伝）。

--- 予言（実行前に登録） ---
genesis/gray_scott_3d.py の mitosis パラメータ域(F=0.0367, k=0.0649、req-d2-grayscott-mitosis
で3D分裂を既に確認済み)に、双安定の遺伝タグ場T（make_seed_initial_multi・step_with_tag、
2D参考実装 T += dt*(0.10*lap(T) + 4.0*V*T*(1-T)*(T-0.5)) の3D版）を追加する。

from-0: 一様(u=1,v=0)+ノイズ種、7個の創始者スポットにそれぞれ独立ランダムな0/1タグ
（初期状態＝第8監査で許可されるfounder-seeds）。分裂の位置・時刻・どちらの娘がどちらの
タグを継承するかは一切命令しない——以降は同じ受動的な拡散+双安定ロック則が続くだけ。

予言: mitosis域では受動的にスポットが分裂増殖する(division_not_seeded)。分裂後も各スポット
のタグは双安定ロックにより0か1にクリーンに留まり(state_inherited、purity>0.6の割合が高い)、
新しいスポットは既存スポットの近傍にのみ出現し総スポット体積が連続的に変化する
(accounting_consistent、複製≠単純コピー)。全3要件が揃えばreached_level=7、division のみ
ならL7-partialと明記する。

決定的対照: タグなし版(T追跡なし、純粋にgray_scott_3d.pyの既存step)では「遺伝」概念自体が
測定不能——state_inherited=Falseとなり、必然的にL7-partial(division_not_seededのみ)になる
ことを確認する。

honest floor: spots ≠ life（Gray-Scottは反応拡散の場現象、Pearson 1993、tier=measured）。
"""

import os
import sys

import numpy as np
from scipy import ndimage

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io, diagnostics as diag  # noqa: E402
from genesis import gray_scott_3d as gs  # noqa: E402
from genesis.solvers import k_grid  # noqa: E402

SHAPE = (64, 64, 64)
STEPS = 6000
N_SEEDS = 7
SEED_RADIUS = 4


def _periodic_delta(c0, c1, shape):
    d = np.abs(np.asarray(c0) - np.asarray(c1))
    d = np.minimum(d, np.asarray(shape) - d)
    return float(np.sqrt((d ** 2).sum()))


def run_tagged(room_id, seed=1, use_tags=True, notes_extra=""):
    p = dict(gs.DEFAULTS)
    rng = np.random.default_rng(seed)
    u, v, T, founder_centers, founder_tags = gs.make_seed_initial_multi(
        SHAPE, rng, n_seeds=N_SEEDS, seed_radius=SEED_RADIUS)
    _, k2 = k_grid(SHAPE)
    mass0 = float(np.mean(u + v))

    snapshot_every = max(1, STEPS // 60)
    snapshots = []
    diverged = False
    for t in range(STEPS):
        if use_tags:
            u, v, T = gs.step_with_tag(u, v, T, gs.DT, p, k2)
        else:
            u, v = gs.step(u, v, gs.DT, p, k2)
        finite_ok = np.all(np.isfinite(u)) and np.all(np.isfinite(v))
        if use_tags:
            finite_ok = finite_ok and np.all(np.isfinite(T))
        if not finite_ok:
            diverged = True
            break
        if t % snapshot_every == 0 or t == STEPS - 1:
            n, labeled, sizes = diag.connected_components(v > 0.1)
            centroids = []
            valid_sizes = []
            for i in range(1, n + 1):
                if sizes[i - 1] >= 5:
                    com = np.array([float(x) for x in ndimage.center_of_mass(labeled == i)])
                    centroids.append(com)
                    valid_sizes.append(float(sizes[i - 1]))
            spots = gs.spot_tag_purity(v, T) if use_tags else []
            snapshots.append({"step": t, "n_spots": len(centroids), "centroids": centroids,
                              "sizes": valid_sizes, "spots": spots,
                              "total_volume": float(sum(valid_sizes))})
    if not snapshots:
        snapshots = [{"step": 0, "n_spots": 0, "centroids": [], "sizes": [], "spots": [],
                      "total_volume": 0.0}]
    mass1 = float(np.mean(u + v))

    # --- division_not_seeded ---
    n_initial = snapshots[0]["n_spots"]
    n_max = max(s["n_spots"] for s in snapshots)
    n_final = snapshots[-1]["n_spots"]
    division_not_seeded = bool(n_max > max(n_initial, 1) and n_final >= 2)

    # --- state_inherited: 最終スナップショットの各スポットのタグがclean(purity>0.6)な割合 ---
    final_spots = snapshots[-1]["spots"]
    n_clean = sum(1 for s in final_spots if s["clean"])
    clean_fraction = (n_clean / len(final_spots)) if final_spots else None
    state_inherited = bool(use_tags and final_spots and clean_fraction is not None
                           and clean_fraction > 0.9)

    # --- accounting_consistent: 複製が「単純コピー」でなく物理的に連続的か ---
    # (a) 局所性: 新規スポットは直前スナップショットの既存スポット近傍にのみ出現するか
    #     (周期距離。無関係な空の場所から突然湧かない=分裂位置を命令していないことの帰結確認)。
    # (b) 連続性: スポット総体積(sizes合計)の相対変化が、分裂ステップで異常な飛躍を示さないか。
    locality_violations = 0
    locality_checks = 0
    for i in range(1, len(snapshots)):
        prev, cur = snapshots[i - 1], snapshots[i]
        if len(cur["centroids"]) <= len(prev["centroids"]) or not prev["centroids"]:
            continue
        for c in cur["centroids"]:
            locality_checks += 1
            dmin = min(_periodic_delta(c, c0, SHAPE) for c0 in prev["centroids"])
            if dmin > 6.0 * SEED_RADIUS:  # 分裂前の伸長+ピンチオフの典型スケールを大きく超える
                locality_violations += 1
    locality_ok = bool(locality_checks == 0 or locality_violations == 0)

    # 体積変化の「通常時」基準は、創始者スポットが安定サイズへ落ち着くまでの初期過渡
    # (warmup、分裂とは無関係な物理的整定であって分裂事故ではない)と、分裂そのものが
    # 起きたステップを除いた「静穏期」(spot数不変の区間)だけから求める——そうしないと
    # 初期過渡の大きな変動に基準がひきずられ、真の分裂ジャンプを検出できなくなる。
    n_spots_series = [s["n_spots"] for s in snapshots]
    vol_series = np.array([s["total_volume"] for s in snapshots])
    rel_changes = np.abs(np.diff(vol_series)) / np.maximum(vol_series[:-1], 1.0)
    warmup = max(3, len(snapshots) // 6)
    quiescent = [i for i in range(warmup, len(snapshots))
                 if n_spots_series[i] == n_spots_series[i - 1]]
    # warmup中の分裂イベントは「創始者スポットが安定サイズへ整定する過程」と本当の自己複製
    # とが混ざり合っており評価不能——warmup後(定常運転下)の分裂のみを会計評価の対象にする。
    division_idx = [i for i in range(warmup, len(snapshots)) if n_spots_series[i] > n_spots_series[i - 1]]
    n_division_in_warmup = sum(1 for i in range(1, warmup) if n_spots_series[i] > n_spots_series[i - 1])
    baseline = rel_changes[[i - 1 for i in quiescent]] if len(quiescent) >= 5 else rel_changes
    if len(baseline) >= 5:
        typical = float(np.median(baseline))
        spread = float(np.std(baseline)) + 1e-9
        division_jumps = rel_changes[[i - 1 for i in division_idx]] if division_idx else np.array([])
        outliers = int(np.sum(division_jumps > typical + 6.0 * spread))
        volume_continuous = bool(outliers == 0)
    else:
        typical, spread, outliers, volume_continuous = None, None, 0, True

    accounting_consistent = bool(locality_ok and volume_continuous)

    reached_L7 = bool(division_not_seeded and state_inherited and accounting_consistent)
    l7_partial = bool(division_not_seeded and not reached_L7)

    detected_keys = ["difference", "localization", "spontaneous_motion", "circulation",
                      "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                      "growth_division_inheritance", "selection_open_ended"]
    detected = {k: False for k in detected_keys}
    detected["difference"] = True
    detected["localization"] = bool(n_final >= 1)
    detected["persistent_individuality"] = bool(n_final >= 1 and not diverged)
    detected["growth_division_inheritance"] = reached_L7

    reached = 2
    if detected["persistent_individuality"]:
        reached = 4
    if division_not_seeded:
        reached = 6  # 分裂のみでは6止まり(partial)、下でreached_L7なら7に昇格
    if reached_L7:
        reached = 7

    n_series = [(s["step"] * gs.DT, s["n_spots"]) for s in snapshots]
    role = "E" if reached_L7 else ("N" if division_not_seeded else "F")
    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {
            "use_tags": use_tags, "F": p["F"], "k": p["k"], "n_seeds": N_SEEDS,
            "founder_tags": founder_tags, "grid": list(SHAPE), "steps": STEPS,
            "n_spots_initial": n_initial, "n_spots_max": n_max, "n_spots_final": n_final,
            "n_series_sampled": n_series[::max(1, len(n_series) // 25)],
            "division_not_seeded": division_not_seeded,
            "final_spot_purities": [s["purity"] for s in final_spots] if final_spots else [],
            "clean_fraction": round(clean_fraction, 4) if clean_fraction is not None else None,
            "state_inherited": state_inherited,
            "locality_checks": locality_checks, "locality_violations": locality_violations,
            "locality_ok": locality_ok,
            "volume_rel_change_typical": round(typical, 6) if typical is not None else None,
            "volume_rel_change_outliers": outliers,
            "volume_continuous": volume_continuous,
            "accounting_warmup_snapshots_excluded": warmup,
            "divisions_in_warmup_excluded_from_accounting": n_division_in_warmup,
            "accounting_consistent": accounting_consistent,
            "l7_gate": {"division_not_seeded": division_not_seeded,
                       "state_inherited": state_inherited,
                       "accounting_consistent": accounting_consistent,
                       "reached_L7": reached_L7},
            "l7_partial": l7_partial,
            "honest_note": ("state_inheritedは「各スポットのタグがclean(0/1にロック)されて"
                            "いるか」の割合(bistable purity>0.6の閾値)で測る——依頼書の定義通り。"
                            "娘スポットが親と同じ系統IDを持つことの明示的なlineage trackingは"
                            "行っていない(親子関係の直接追跡でなく、事後のタグ純度で代理測定)。"
                            "accounting_consistentは(a)新規スポットが直前の既存スポット近傍"
                            "(周期距離)にのみ出現するか(空の場所からの湧き出しでないか)、"
                            "(b)スポット総体積の相対変化に分裂時の異常な飛躍がないか、の2条件。"
                            "(b)の「通常時」基準は創始者スポットが安定サイズへ整定するまでの"
                            "初期過渡(warmup、accounting_warmup_snapshots_excluded)を除外して"
                            "求める——この間の体積変化は自己複製でなく初期条件からの物理的な"
                            "整定であり、除外しないと真の分裂ジャンプの基準がゆがむ。同じ理由で"
                            "warmup中に起きた分裂(divisions_in_warmup_excluded_from_accounting)"
                            "はaccounting評価の対象外とする(整定と複製が混ざり評価不能なため)。"),
        },
        "purity": {"per_object_labels": False, "external_optimum": False, "role": role},
    }
    integrity = io.integrity_block(
        conservation_drift=mass1 - mass0, resolutions_result={"%dx%dx%d" % SHAPE: n_max},
        seed_success={str(seed): reached}, nan_or_clip=diverged)
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=False,  # タグは創始者スポットの初期状態のみ、分裂/継承先は指定しない
        gate_encodes_conclusion_causality=False,
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "no-tag control (same F,k, T追跡なし)",
                        "result": "see companion room -- state_inherited測定不能=L7-partialのはず"}])
    checksum = io.checksum_of([u, v] + ([T] if use_tags else []))
    genesis_yaml = {
        "equations": "Gray-Scott (du/dt=Du*lap(u)-u*v^2+F(1-u), dv/dt=Dv*lap(v)+u*v^2-(F+k)v) "
                     "+ 双安定遺伝タグ場 dT/dt=Dv*lap(T)+4*v*T*(1-T)*(T-0.5) (vでゲート)"
                     if use_tags else
                     "Gray-Scott (du/dt=Du*lap(u)-u*v^2+F(1-u), dv/dt=Dv*lap(v)+u*v^2-(F+k)v), タグなし対照",
        "solver": "pseudo-spectral, semi-implicit diffusion, explicit reaction",
        "dt": gs.DT, "dx": 1.0, "grid": list(SHAPE), "boundary": "periodic",
        "params": p, "seed": seed, "seeds": [seed], "commit": None, "checksum": checksum,
    }
    notes = ("確認B [%s] F=%.4f k=%.4f, steps=%d, n_seeds=%d, founder_tags=%s。"
             "n_spots: %d->%d(max=%d)。division_not_seeded=%s, clean_fraction=%s, "
             "state_inherited=%s, locality_ok=%s(%d/%d), volume_continuous=%s(outliers=%d), "
             "accounting_consistent=%s -> reached_L7=%s(l7_partial=%s)。%s"
             % (notes_extra, p["F"], p["k"], STEPS, N_SEEDS, founder_tags, n_initial, n_final,
                n_max, division_not_seeded, clean_fraction, state_inherited, locality_ok,
                locality_violations, locality_checks, volume_continuous, outliers,
                accounting_consistent, reached_L7, l7_partial, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  n_spots %d->%d(max=%d) division=%s state_inherited=%s "
          "accounting=%s reached_L7=%s(partial=%s)"
          % (run_dir, n_initial, n_final, n_max, division_not_seeded, state_inherited,
             accounting_consistent, reached_L7, l7_partial))
    return snapshots, report


def main():
    print("=== 確認B: 完全なL7を3D tagged Gray-Scottで（分裂＋遺伝） ===")
    print("予言: mitosis域で受動分裂 AND タグがclean(purity>0.6)に留まる AND 体積会計が連続 "
          "-> reached_L7。タグなし対照はL7-partial止まり。")

    print("\n[タグあり: mitosis域] F=0.0367 k=0.0649, 7創始者スポット")
    run_tagged("l7-tagged-mitosis-seed0001", seed=1, use_tags=True,
              notes_extra="タグあり: 分裂+遺伝+会計の全3要件を期待。")

    print("\n[決定的対照: タグなし] 同一F,k、Tを追跡しない")
    run_tagged("l7-notag-control-seed0001", seed=1, use_tags=False,
              notes_extra="決定的対照(タグなし): 分裂は起きるが遺伝という概念自体が測定不能"
                          "->L7-partial(division_not_seededのみ)になることを期待。")

    print("=== 確認B done ===")


if __name__ == "__main__":
    main()
