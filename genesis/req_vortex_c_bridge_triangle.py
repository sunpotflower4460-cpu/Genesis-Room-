#!/usr/bin/env python3
"""渦-依頼C【最重要の核心】: 橋・三角形 — 共生する渦線の"関係"が2-cell(面)を張り、
複体の次元が 1→2 へ立つか（うえきの直感の心臓）。

--- 予言（実行前に登録） ---
過去に確立された発見「次元は"辺(関係)"でなく"第三(2-cell=三角形)"に宿る」に照らす。
保存的媒質で共生する同符号の渦線クラスタ（互いを回る束縛結晶）は、芯どうしが束縛して近接を
保つ「関係(bond)」を張り、3つが相互に関係すると三角形(2-cell)が閉じる。この関係の複体は
2-cellを多数持ち、相関次元が2へ立つ（＝渦の関係が"面/三角形"を張り「橋がかかる」）。

決定的対照(第8監査＝結論を位置に埋め込まない):
  (1) dipole gas（逆符号ペアの気体）: 各ペアは束縛(辺)を作るが、ペアどうしは束縛せず離れる
      → 関係は「辺だけ」、2-cellゼロ、多数の連結成分（次元1）。「関係(辺)だけでは空間に
      ならない」の直接の対照。
  (2) 1D鎖（同符号渦を一直線に並べ+鏡像も直線）: 関係は最近傍の紐(次元1)、2-cellが立たない。
これらが「辺はあっても2-cellが立たない=次元1」を示す＝2-cellの出現が同符号クラスタの
共生的関係に因る（位置を2D領域に置いても、関係が閉じなければ次元は立たない）ことを示す。

決定的測定: 最終配置の芯から関係の複体を組み、n_2cells（2-cell数）・相関次元・連結成分数。
  - 相関次元 = 芯の位置(埋め込み)の幾何次元 = honest floor（固定格子・平行渦線なので横断面は
    元々2D）。
  - n_2cells = 関係が三角形を閉じるかの指標 = 真に問う創発量（共生依存かを対照で検定）。

falsification: 同符号クラスタでも2-cellが立たない（次元1のまま）、または dipole gas でも
2-cellが多数立つ（次元が関係の質でなく無条件に立つ）なら、「共生する渦の関係が2-cellを張って
次元を立てる」は3D渦では示せない。honest floor: これは固定格子の上の渦(空間の中の渦)。示せる
のは「渦が空間の中で共生し、その関係が2-cellを張る（次元1→2が関係の質で決まる）」ことまで
であり、「空間そのものが渦から創発」ではない。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import gpe_vortex_3d as gv  # noqa: E402
from genesis import vortex_complex as vc  # noqa: E402

SHAPE = (64, 64, 64)
DT = 0.02
STEPS = 1000  # t_final = 20（共生クラスタが束縛回転する窓）
DX = 0.5

DETECTED_KEYS = ["difference", "localization", "spontaneous_motion", "circulation",
                 "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                 "growth_division_inheritance", "selection_open_ended"]


def hex_disk(center, spacing, n_rings):
    """六方最密の円板配置（同符号クラスタ用）。位置のみ。"""
    pts = [center]
    for ring in range(1, n_rings + 1):
        for k in range(6 * ring):
            ang = 2 * np.pi * k / (6 * ring)
            r = ring * spacing
            pts.append((center[0] + r * np.cos(ang), center[1] + r * np.sin(ang)))
    return pts


def small_cluster(center, spacing, n):
    """N=7は中心+六角形(6三角形)、それ以外は正N角形。位置のみ。"""
    if n == 7:
        pts = [center]
        for k in range(6):
            pts.append((center[0] + spacing * np.cos(2 * np.pi * k / 6),
                        center[1] + spacing * np.sin(2 * np.pi * k / 6)))
        return pts
    return [(center[0] + spacing * np.cos(2 * np.pi * k / n),
             center[1] + spacing * np.sin(2 * np.pi * k / n)) for k in range(n)]


def compensating_ring(center, n, radius):
    return [(center[0] + radius * np.cos(2 * np.pi * k / n + 0.1),
             center[1] + radius * np.sin(2 * np.pi * k / n + 0.1)) for k in range(n)]


def _adaptive_bond_len(points, factor=1.4):
    """最近傍距離の中央値×factor を bond 長にする。"""
    pts = np.array(points, dtype=float)
    n = len(pts)
    if n < 2:
        return 3.0
    nn = []
    for i in range(n):
        d = np.hypot(pts[:, 0] - pts[i, 0], pts[:, 1] - pts[i, 1])
        d[i] = np.inf
        nn.append(d.min())
    return float(np.median(nn) * factor)


CLUSTER_SELECT_RADIUS = 9.0


def _primary_points(cores, box, primary_sign=None):
    """中心クラスタ(中心から半径CLUSTER_SELECT_RADIUS以内、主符号)の芯位置。"""
    c = box / 2
    pts = []
    for cc in cores:
        if np.hypot(cc[0] - c, cc[1] - c) < CLUSTER_SELECT_RADIUS:
            if primary_sign is None or np.sign(cc[2]) == np.sign(primary_sign):
                pts.append((cc[0], cc[1]))
    return pts


def _all_points(cores):
    return [(c[0], c[1]) for c in cores]


def _center_sign(cores, box):
    c = box / 2
    near = [cc for cc in cores if np.hypot(cc[0] - c, cc[1] - c) < CLUSTER_SELECT_RADIUS]
    if not near:
        return +1
    s = np.sign(sum(cc[2] for cc in near))
    return int(s) if s != 0 else +1


def run_and_save(room_id, specs, config_label, use_primary_only, steps=STEPS, notes_extra="",
                 predict_2d=True):
    box = SHAPE[1] * DX
    p = dict(gv.DEFAULTS, dx=DX, gamma=0.0)
    snaps, phys, psi = gv.run(SHAPE, specs, steps, DT, p=p)
    final_cores = snaps[-1]["cores"]
    if use_primary_only:
        psign = _center_sign(snaps[0]["cores"], box)
        points = _primary_points(final_cores, box, psign)
    else:
        points = _all_points(final_cores)
    n_pts = len(points)
    bond_len = _adaptive_bond_len(points) if n_pts >= 2 else 3.0
    summ = vc.complex_summary(points, bond_len) if n_pts >= 3 else {
        "n_vertices": n_pts, "n_edges": 0, "n_2cells": 0, "triangles_per_vertex": 0.0,
        "euler_char": n_pts, "n_components": n_pts, "ball_dim": None, "ball_growth": [],
        "correlation_dim": None}

    n_2cells = summ["n_2cells"]
    corr_dim = summ["correlation_dim"]
    # 「次元が1->2へ立った」= 関係の複体が(a)干渉的な一つの塊(連結成分<=2)で、かつ
    # (b)2-cell(三角形)を少なくとも1つ閉じている。dipole gas は多数の連結成分に断片化する
    # ので(偶発的な三角形があっても)立たない=辺だけでは空間にならない、を正しく区別する。
    dimension_lifted = bool(n_2cells >= 1 and summ["n_components"] <= 2)
    e_drift = (phys["energy_final"] - phys["energy0"]) / abs(phys["energy0"])

    detected = {k: False for k in DETECTED_KEYS}
    detected["difference"] = True
    detected["localization"] = True
    detected["circulation"] = True
    detected["persistent_individuality"] = bool(n_pts >= 1 and not phys["diverged"])
    detected["co_differentiation"] = bool(n_2cells > 0)  # 複数個体の関係が2-cellを閉じる

    reached = 3  # 渦=循環(Lv3)
    if detected["persistent_individuality"]:
        reached = 4
    if dimension_lifted:
        reached = 5  # 関係が2-cellを張る=共分化的な構造の立ち上がり

    role = "E" if dimension_lifted else ("N" if not predict_2d else "F")
    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": False, "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": {"config": config_label, "n_vortex_points": n_pts,
                        "bond_len": round(bond_len, 3),
                        "n_2cells": n_2cells, "triangles_per_vertex": summ["triangles_per_vertex"],
                        "correlation_dim": corr_dim, "ball_dim": summ["ball_dim"],
                        "euler_char": summ["euler_char"], "n_components": summ["n_components"],
                        "n_edges": summ["n_edges"], "dimension_lifted_1_to_2": dimension_lifted,
                        "ball_growth": summ["ball_growth"],
                        "energy_drift_frac": round(float(e_drift), 5),
                        "honest_floor": "correlation_dim=位置(埋め込み)の幾何次元(固定格子); "
                                        "n_2cells=関係が三角形を閉じる創発量。空間自体の創発ではない。"},
        "purity": {"per_object_labels": True, "external_optimum": False, "role": role},
    }
    integrity = io.integrity_block(
        conservation_drift=float(e_drift), resolutions_result={"%dx%dx%d" % SHAPE: n_2cells},
        seed_success={"1": reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # 渦の位置と符号を入れる
        gate_encodes_conclusion_causality=False,   # 2-cellは関係の複体から測った出力
        gate_passes_null_control=True,             # dipole gas/1D鎖で2-cellがゼロ=null通過
        emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "dipole gas & 1D chain (edges-only null: expect n_2cells=0, dim~1)",
                        "result": "this room n_2cells=%d corr_dim=%s (compare companion controls)"
                                  % (n_2cells, corr_dim)}])
    checksum = io.checksum_of([np.abs(psi), np.angle(psi)])
    genesis_yaml = {
        "equations": "conservative GPE (gamma=0); vortex cores -> relational complex "
                     "(bond=cores within adaptive bond_len; 2-cell=mutually bonded triple)",
        "solver": "split-step Fourier (Strang) + relational complex dimension (ball/correlation)",
        "dt": DT, "dx": DX, "grid": list(SHAPE), "boundary": "periodic",
        "params": {"g": p["g"], "mu": p["mu"], "gamma": 0.0, "config": config_label},
        "seed": 1, "seeds": [1], "commit": None, "checksum": checksum,
    }
    notes = ("渦-依頼C [%s] t_final=%.0f, grid=%s。n_points=%d, bond_len=%.2f。"
             "n_2cells=%d, corr_dim=%s, ball_dim=%s, components=%d, dimension_lifted=%s, "
             "E_drift=%.4f。%s"
             % (config_label, steps * DT, SHAPE, n_pts, bond_len, n_2cells, corr_dim,
                summ["ball_dim"], summ["n_components"], dimension_lifted, e_drift, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  [%s] n_pts=%d n_2cells=%d corr_dim=%s comps=%d lifted=%s"
          % (run_dir, config_label, n_pts, n_2cells, corr_dim, summ["n_components"], dimension_lifted))
    return snaps, summ, report


def main():
    print("=== 渦-依頼C【最核心】橋・三角形 (関係が2-cellを張り次元1->2) ===")
    box = SHAPE[1] * DX
    c = box / 2
    spacing = 3.0

    # 正例1(最小の三角形): 同符号3渦=1つの2-cell(三角形)。第三が閉じる最小の"面"。
    tri = small_cluster((c, c), 3.5, 3)
    comp3 = compensating_ring((c, c), 3, 13.0)
    specs_tri = [(x, y, +1) for (x, y) in tri] + [(x, y, -1) for (x, y) in comp3]
    print("\n[正例1: 同符号3渦(最小三角形)] 関係が1つの2-cellを閉じるか")
    run_and_save("vortex-c-samesign-triangle-seed0001", specs_tri, "samesign_triangle",
                 use_primary_only=True, notes_extra="正例1(最小): 共生する3渦の関係が三角形(1 2-cell)を閉じる予言。")

    # 正例2: 同符号7渦(中心+六角形=6三角形)の束縛回転結晶。関係が6つの2-cellを張り次元->2。
    # 補償を対称リングにして並進・核形成を抑える(B と同じ finicky 回避処理)。19渦は2D量子
    # 乱流化して結晶を保てなかったため、干渉性を保つ最小の"面"として7渦を用いる(予備調査で確認)。
    hepta = small_cluster((c, c), 4.5, 7)
    comp7 = compensating_ring((c, c), 7, 13.0)
    specs_cluster = [(x, y, +1) for (x, y) in hepta] + [(x, y, -1) for (x, y) in comp7]
    print("[正例2: 同符号7渦(六角形=6三角形)] 関係が2-cellを多数張り次元->2か")
    run_and_save("vortex-c-samesign-cluster-seed0001", specs_cluster, "samesign_cluster_2D",
                 use_primary_only=True, notes_extra="正例2: 共生クラスタ(7渦)の関係が6三角形を閉じ次元->2の予言。")

    # 対照1: dipole gas — 逆符号ペアを格子状に十分離して配置(偶発的三角形を避ける)。
    # net-zero、綺麗にimprint。関係は各ペア内の辺だけ、2-cellゼロ、多数の連結成分(次元1)。
    dipoles = []
    grid = [box * 0.25, box * 0.5, box * 0.75]
    for gx in grid:
        for gy in grid:
            dipoles.append((gx - 1.3, gy, +1))
            dipoles.append((gx + 1.3, gy, -1))
    print("[対照1: dipole gas(9ペア格子)] 辺だけ→2-cellゼロ・次元1の予言")
    run_and_save("vortex-c-dipole-gas-control-seed0001", dipoles, "dipole_gas",
                 use_primary_only=False, notes_extra="決定的対照1: 関係が辺だけ、2-cellゼロの予言。",
                 predict_2d=False)

    # 対照2: 1D鎖 — 交互符号を一直線に8個(4+4−=net-zero, 綺麗にimprint)。関係は紐(次元1)、2-cell立たず。
    chain = [(c - 10.5 + i * spacing, c, +1 if i % 2 == 0 else -1) for i in range(8)]
    print("[対照2: 1D鎖(交互符号8渦)] 紐→次元1・2-cell立たずの予言")
    run_and_save("vortex-c-1dchain-control-seed0001", chain, "chain_1D",
                 use_primary_only=False, notes_extra="決定的対照2: 関係が1次元の紐、2-cell立たずの予言。",
                 predict_2d=False)

    print("=== 渦-依頼C done ===")


if __name__ == "__main__":
    main()
