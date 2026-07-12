#!/usr/bin/env python3
"""渦-依頼B【核心】: 共生しながら"登る" — 束縛状態を作り、第三が結合、各々巻きを保つ。

--- 予言（実行前に登録） ---
保存的媒質(γ=0)で、同符号の渦線が「互いを回る軌道」(2Dで確立)を3Dで示し、
  (i) 2つの同符号渦 → 束縛"分子"として co-rotate（分子として持続、離散せず）。
  (ii) 3つ → 束縛三角形として co-rotate（2→3の"登り"）。第三が結合して三体束縛。
  (iii) 4つ → 束縛四辺形。
  いずれも各渦が自分の整数巻き(+1)を保つ＝同一化(単一の+N芯へ併合)せず、N個の識別可能な
  芯が持続する。総巻き数を厳密保存。「同一化せず共存しながら複雑化」。
周期箱は net-zero を要求するため、主クラスタ(+1×N、左)の反対側遠方に鏡像クラスタ(−1×N、右)を
置き、主クラスタ(x<box/2)だけを追跡する。

決定的測定: 主クラスタの芯数 n_primary が N を維持（併合しない＝アイデンティティ維持）、
クラスタが束縛（芯間の最大距離が有界＝離散しない分子）、co-rotate（累積回転>0）、巻き保存。

決定的対照(第8監査): 同一配置で γ を 0→有限（散逸）にすると、束縛分子は崩れる/併合する
（多重アイデンティティの持続には保存的媒質が要る）ことを対照として示す。

falsification: 保存的でも同符号クラスタが即座に単一+N芯へ併合する（識別可能な複数が持続
しない）、または第三が常に射出されて束縛三角形が作れない（登れない）なら、「共生しながら
登る」は3D渦線で崩れる。honest floor: 固定格子上の渦。示せるのは「渦が空間の中で束縛分子を
作り、識別性を保って積み上がる」ことまで。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import gpe_vortex_3d as gv  # noqa: E402

SHAPE = (64, 64, 64)
DT = 0.02
STEPS = 1500  # t_final = 30
DX = 0.5

DETECTED_KEYS = ["difference", "localization", "spontaneous_motion", "circulation",
                 "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                 "growth_division_inheritance", "selection_open_ended"]


def polygon(center, radius, n, phase0=0.0):
    if n == 1:
        return [(center[0], center[1])]
    return [(center[0] + radius * np.cos(phase0 + 2 * np.pi * i / n),
             center[1] + radius * np.sin(phase0 + 2 * np.pi * i / n)) for i in range(n)]


PRIM_RADIUS_SELECT = 7.0  # 中心クラスタとみなす半径
COMP_RING_RADIUS = 11.0   # 補償リング半径


def cluster_with_mirror(box, n, cluster_radius=2.4):
    """中心の同符号クラスタ(+1×N) + 対称な補償リング(−1×N、大半径)。
    補償を対称に配置して正味の双極子モーメントを消す(=クラスタが並進せず、芯どうしが箱内で
    衝突・核形成しない)ことが、tight な同符号クラスタを綺麗に保つ鍵——周期箱の finicky さ
    (Claudeの3D場sim が濁った件)を回避する処理。位置と符号のみを入れる。"""
    c = (box * 0.5, box * 0.5)
    specs = []
    for (x, y) in polygon(c, cluster_radius, n, phase0=0.0):
        specs.append((x, y, +1))
    for (x, y) in polygon(c, COMP_RING_RADIUS, n, phase0=np.pi / max(n, 1)):
        specs.append((x, y, -1))
    return specs


def _primary_sign(snaps, box):
    """初期フレームで中心近傍(r<PRIM_RADIUS_SELECT)に居る芯の多数派符号を主クラスタの符号とする。"""
    c = box * 0.5
    first = snaps[0]["cores"]
    near = [cc for cc in first if np.hypot(cc[0] - c, cc[1] - c) < PRIM_RADIUS_SELECT]
    if not near:
        return +1
    s = np.sign(sum(cc[2] for cc in near))
    return int(s) if s != 0 else +1


def _primary_stats(cores, box, primary_sign=None):
    """中心クラスタの芯(中心から半径PRIM_RADIUS_SELECT以内かつ主符号): 個数, 芯間最大距離, 重心。"""
    c = box * 0.5
    if primary_sign is None:
        prim = [cc for cc in cores if np.hypot(cc[0] - c, cc[1] - c) < PRIM_RADIUS_SELECT]
    else:
        prim = [cc for cc in cores if np.hypot(cc[0] - c, cc[1] - c) < PRIM_RADIUS_SELECT
                and np.sign(cc[2]) == np.sign(primary_sign)]
    n = len(prim)
    if n == 0:
        return {"n": 0, "max_dist": None, "centroid": None, "cores": []}
    maxd = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            d = np.hypot(prim[i][0] - prim[j][0], prim[i][1] - prim[j][1])
            maxd = max(maxd, d)
    cen = (float(np.mean([c[0] for c in prim])), float(np.mean([c[1] for c in prim])))
    return {"n": n, "max_dist": round(float(maxd), 3), "centroid": cen, "cores": prim}


def _cum_rotation(snaps, box, primary_sign):
    cum = 0.0
    prev = _primary_stats(snaps[0]["cores"], box, primary_sign)
    for s in snaps[1:]:
        cur = _primary_stats(s["cores"], box, primary_sign)
        if prev["n"] >= 2 and cur["n"] >= 2 and prev["centroid"] and cur["centroid"]:
            cen = prev["centroid"]
            deltas = []
            for (x, y, q, *_r) in prev["cores"]:
                nn = min(cur["cores"], key=lambda cc: np.hypot(cc[0] - x, cc[1] - y))
                a0 = np.arctan2(y - cen[1], x - cen[0])
                a1 = np.arctan2(nn[1] - cen[1], nn[0] - cen[0])
                deltas.append(np.degrees((a1 - a0 + np.pi) % (2 * np.pi) - np.pi))
            if deltas:
                cum += float(np.mean(deltas))
        prev = cur
    return round(cum, 2)


def run_and_save(room_id, n, gamma, steps=STEPS, cluster_radius=2.2, notes_extra=""):
    box = SHAPE[1] * DX
    specs = cluster_with_mirror(box, n, cluster_radius=cluster_radius)
    p = dict(gv.DEFAULTS, dx=DX, gamma=gamma)
    snaps, phys, psi = gv.run(SHAPE, specs, steps, DT, p=p)

    primary_sign = _primary_sign(snaps, box)
    prim0 = _primary_stats(snaps[0]["cores"], box, primary_sign)
    primf = _primary_stats(snaps[-1]["cores"], box, primary_sign)
    n_prim_series = [_primary_stats(s["cores"], box, primary_sign)["n"] for s in snaps]
    maxd_series = [_primary_stats(s["cores"], box, primary_sign)["max_dist"] for s in snaps]
    identity_kept = bool(all(x == n for x in n_prim_series))  # 常にN個の識別可能な芯(符号で追跡)
    finite_maxd = [d for d in maxd_series if d is not None]
    # 束縛=同符号クラスタの芯間最大距離が有界(初期クラスタ径の数倍以内、離散しない分子)
    bounded = bool(finite_maxd and max(finite_maxd) < 4.0 * cluster_radius + 3.0)
    cum_rot = _cum_rotation(snaps, box, primary_sign)
    net_windings = [s["net_winding"] for s in snaps]
    winding_conserved = bool(all(w == net_windings[0] for w in net_windings))
    e_drift = (phys["energy_final"] - phys["energy0"]) / abs(phys["energy0"])
    molecule_persists = bool(identity_kept and bounded and not phys["diverged"])

    detected = {k: False for k in DETECTED_KEYS}
    detected["difference"] = True
    detected["localization"] = True
    detected["circulation"] = True
    detected["spontaneous_motion"] = bool(abs(cum_rot) > 5.0)
    detected["persistent_individuality"] = molecule_persists and winding_conserved
    # co_differentiation(Lv5)的: 複数の識別可能な個体が関係を保ち共存(併合せず)
    detected["co_differentiation"] = bool(n >= 3 and molecule_persists)

    reached = 2
    if detected["circulation"] or detected["spontaneous_motion"]:
        reached = 3
    if detected["persistent_individuality"]:
        reached = 4
    if detected["co_differentiation"]:
        reached = 5

    role = "E" if molecule_persists else ("N" if gamma > 0 else "F")
    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": False, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"n_same_sign": n, "gamma": gamma,
                        "n_primary_initial": prim0["n"], "n_primary_final": primf["n"],
                        "identity_kept_N_throughout": identity_kept,
                        "bound_molecule_persists": molecule_persists,
                        "max_intra_distance_max": round(max(finite_maxd), 3) if finite_maxd else None,
                        "cluster_radius0": cluster_radius,
                        "cum_rotation_deg": cum_rot,
                        "winding_conserved": winding_conserved,
                        "energy_drift_frac": round(float(e_drift), 5),
                        "n_primary_series": list(zip([s["t"] for s in snaps], n_prim_series))[::max(1, len(snaps) // 12)],
                        "max_dist_series": list(zip([s["t"] for s in snaps], maxd_series))[::max(1, len(snaps) // 12)]},
        "purity": {"per_object_labels": True, "external_optimum": False, "role": role},
    }
    integrity = io.integrity_block(
        conservation_drift=float(e_drift), resolutions_result={"%dx%dx%d" % SHAPE: primf["n"]},
        seed_success={"1": reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # 渦の位置と符号を入れる
        gate_encodes_conclusion_causality=False,  # 束縛/併合は測定した出力
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "dissipative (gamma>0) companion breaks the molecule",
                        "result": "molecule_persists=%s at gamma=%.2f" % (molecule_persists, gamma)}])
    checksum = io.checksum_of([np.abs(psi), np.angle(psi)])
    genesis_yaml = {
        "equations": "(i - gamma) dpsi/dt = [-1/2 lap + g|psi|^2 - mu] psi (damped GPE); "
                     "N same-sign vortex lines (primary) + N mirror (net-zero periodic box)",
        "solver": "split-step Fourier (Strang), periodic BC",
        "dt": DT, "dx": DX, "grid": list(SHAPE), "boundary": "periodic",
        "params": {"g": p["g"], "mu": p["mu"], "gamma": gamma, "n_same_sign": n},
        "seed": 1, "seeds": [1], "commit": None, "checksum": checksum,
    }
    notes = ("渦-依頼B [N=%d, gamma=%.2f] t_final=%.0f, grid=%s。identity_kept=%s, "
             "molecule_persists=%s, cum_rot=%.1fdeg, winding_conserved=%s, E_drift=%.4f。"
             "n_primary %d->%d。%s"
             % (n, gamma, steps * DT, SHAPE, identity_kept, molecule_persists, cum_rot,
                winding_conserved, e_drift, prim0["n"], primf["n"], notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  [N=%d g=%.2f] identity=%s molecule=%s cum_rot=%.0f Edrift=%.4f nprim %d->%d"
          % (run_dir, n, gamma, identity_kept, molecule_persists, cum_rot, e_drift,
             prim0["n"], primf["n"]))
    return snaps, phys, report


def main():
    print("=== 渦-依頼B【核心】共生しながら登る (束縛分子+第三が結合, 巻き保存) ===")

    print("\n[二体分子・保存的] N=2 同符号 → 束縛co-rotate")
    run_and_save("vortex-b-diatomic-conservative-seed0001", 2, 0.0,
                 notes_extra="二体束縛分子(共生ペア)。")
    print("[三体分子・保存的] N=3 → 束縛三角形(2->3の登り)")
    run_and_save("vortex-b-triatomic-conservative-seed0001", 3, 0.0,
                 notes_extra="三体束縛(第三が結合)。2->3の登り。")
    print("[四体分子・保存的] N=4 → 束縛四辺形")
    run_and_save("vortex-b-tetratomic-conservative-seed0001", 4, 0.0,
                 notes_extra="四体束縛。さらに積む。")
    print("\n[決定的対照・散逸] N=3, gamma=0.3 → 分子が崩れる/併合するか")
    run_and_save("vortex-b-triatomic-dissipative-control-seed0001", 3, 0.30,
                 notes_extra="決定的対照: 散逸下で多重アイデンティティ分子が持続するか。")

    print("=== 渦-依頼B done ===")


if __name__ == "__main__":
    main()
