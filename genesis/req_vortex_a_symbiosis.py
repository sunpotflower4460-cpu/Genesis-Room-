#!/usr/bin/env python3
"""渦-依頼A: クリーンな3D渦共生（媒質依存を綺麗に）。

--- 予言（実行前に登録） ---
2Dで確立した「保存的=共生 vs 散逸=対消滅」を3Dの渦線で再現する。
  保存的媒質(γ=0):
    (i) 逆符号ペア(+1/−1、周期箱で唯一許される net-zero) → 両方生存し並進(共存)。
    (ii) net-zero 4渦(交互符号の正方形) → 4芯すべて生存し束縛回転結晶として回る(共生)。
    いずれも総巻き数=0を厳密保存、GPEエネルギー保存。
  散逸媒質(γ>0):
    逆符号ペア → 芯が埋まって対消滅(n_cores→0、E緩和、min密度が背景へ回復)。

決定的対照(第8監査): 同一の初期渦配置のまま γ を 0→有限 に変えるだけで、
共存 → 対消滅 に切り替わる = 運命を決めるのは初期条件でなく「媒質(保存的か散逸的か)」。

falsification: γ=0でも対消滅する、またはγ>0でも逆符号ペアが共存し続ける(媒質非依存)なら、
「共生の媒質依存性」は3Dで崩れる。honest floor: これは固定格子上の渦(空間の中の渦)であり、
「渦から空間が創発」ではない。ここで示せるのは渦が媒質の性質に従って共存/対消滅すること。
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from common import io  # noqa: E402
from genesis import gpe_vortex_3d as gv  # noqa: E402

SHAPE = (64, 64, 64)
DT = 0.02
STEPS = 1200  # t_final = 24
DX = 0.5

DETECTED_KEYS = ["difference", "localization", "spontaneous_motion", "circulation",
                 "persistent_individuality", "co_differentiation", "self_maintaining_closure",
                 "growth_division_inheritance", "selection_open_ended"]


def _centroid(cores):
    return (float(np.mean([c[0] for c in cores])), float(np.mean([c[1] for c in cores])))


def _core_motion(snaps):
    """芯の運動を要約: 重心の並進距離と、フレーム間で芯を最近傍マッチングして累積した
    符号付き回転角(deg)。対称な回転結晶では mean|angle| は不変で回転を見落とすため、
    連続フレーム間の実際の角変位を積算する。"""
    first = snaps[0]["cores"]
    last = snaps[-1]["cores"]
    if not first or not last:
        return {"translated": None, "cum_rotation_deg": None}
    cx0, cy0 = _centroid(first)
    cxf, cyf = _centroid(last)
    translated = float(np.hypot(cxf - cx0, cyf - cy0))

    cum_rot = 0.0
    for a, b in zip(snaps[:-1], snaps[1:]):
        ca, cb = a["cores"], b["cores"]
        if not ca or not cb:
            continue
        cen = _centroid(ca)
        deltas = []
        for (x, y, q, *_r) in ca:
            # 最近傍の次フレーム芯へマッチ(同符号優先)
            cands = [cc for cc in cb if np.sign(cc[2]) == np.sign(q)] or list(cb)
            nn = min(cands, key=lambda cc: np.hypot(cc[0] - x, cc[1] - y))
            ang0 = np.arctan2(y - cen[1], x - cen[0])
            ang1 = np.arctan2(nn[1] - cen[1], nn[0] - cen[0])
            d = np.degrees((ang1 - ang0 + np.pi) % (2 * np.pi) - np.pi)
            deltas.append(d)
        if deltas:
            cum_rot += float(np.mean(deltas))
    return {"translated": round(translated, 3), "cum_rotation_deg": round(cum_rot, 2)}


def run_and_save(room_id, specs, gamma, config_label, steps=STEPS, notes_extra=""):
    p = dict(gv.DEFAULTS, dx=DX, gamma=gamma)
    snaps, phys, psi = gv.run(SHAPE, specs, steps, DT, p=p)
    n0 = snaps[0]["n_cores"]
    nf = snaps[-1]["n_cores"]
    net_windings = [s["net_winding"] for s in snaps]
    winding_conserved = bool(all(w == net_windings[0] for w in net_windings))
    coexisted = bool(nf >= n0 and nf > 0 and not phys["diverged"])
    annihilated = bool(nf == 0 and n0 > 0)
    e_drift_frac = (phys["energy_final"] - phys["energy0"]) / abs(phys["energy0"])
    motion = _core_motion(snaps)

    detected = {k: False for k in DETECTED_KEYS}
    detected["difference"] = bool(n0 > 0)
    detected["localization"] = bool(n0 > 0)  # 位相巻き=局在トポロジー欠陥
    detected["spontaneous_motion"] = bool(motion["translated"] and motion["translated"] > 0.5)
    detected["circulation"] = bool(n0 > 0)  # 渦=循環そのもの
    detected["persistent_individuality"] = coexisted and winding_conserved

    reached = 0
    if detected["difference"]:
        reached = 1
    if detected["localization"]:
        reached = 2
    if detected["spontaneous_motion"] or detected["circulation"]:
        reached = 3
    if detected["persistent_individuality"]:
        reached = 4

    role = "E" if (coexisted or annihilated) else "F"
    report = {
        "reached_level": reached, "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": False,  # 渦の位置と符号(差の種)は初期条件に置く
        "level_detected_by_measurement": True, "detected": detected,
        "measured_by": {"config": config_label, "gamma": gamma,
                        "n_cores_initial": n0, "n_cores_final": nf,
                        "coexisted": coexisted, "annihilated": annihilated,
                        "net_winding_series": net_windings[::max(1, len(net_windings) // 12)],
                        "winding_conserved": winding_conserved,
                        "energy0": round(phys["energy0"], 3),
                        "energy_final": round(phys["energy_final"], 3),
                        "energy_drift_frac": round(float(e_drift_frac), 5),
                        "norm0": round(phys["norm0"], 3),
                        "norm_final": round(phys["norm_final"], 3),
                        "core_motion": motion,
                        "n_cores_series": [(s["t"], s["n_cores"]) for s in snaps][::max(1, len(snaps) // 15)],
                        "min_density_final": round(snaps[-1]["min_density"], 4)},
        "purity": {"per_object_labels": True, "external_optimum": False, "role": role},
    }
    integrity = io.integrity_block(
        conservation_drift=float(e_drift_frac),
        resolutions_result={"%dx%dx%d" % SHAPE: nf},
        seed_success={"1": reached}, nan_or_clip=phys["diverged"])
    input_vs_output = io.input_output_selfcheck(
        target_encoded_in_initial_condition=True,  # 渦の位置と符号は入れる
        gate_encodes_conclusion_causality=False,  # 共存/対消滅は測定した出力、ゲートに埋めていない
        gate_passes_null_control=False, emerged_quantity_is_algebraic_restatement=False,
        control_runs=[{"name": "medium: gamma 0->finite flips fate (decisive control)",
                        "result": "coexisted=%s annihilated=%s at gamma=%.2f" % (coexisted, annihilated, gamma)}])
    checksum = io.checksum_of([np.abs(psi), np.angle(psi)])
    genesis_yaml = {
        "equations": "(i - gamma) dpsi/dt = [-1/2 lap + g|psi|^2 - mu] psi  (damped GPE); "
                     "gamma=0 conservative(unitary), gamma>0 dissipative(relax to ground state)",
        "solver": "split-step Fourier (Strang), periodic BC",
        "dt": DT, "dx": DX, "grid": list(SHAPE), "boundary": "periodic",
        "params": {"g": p["g"], "mu": p["mu"], "gamma": gamma}, "seed": 1, "seeds": [1],
        "commit": None, "checksum": checksum,
    }
    notes = ("渦-依頼A [%s] gamma=%.2f, t_final=%.0f, grid=%s。n_cores %d->%d, coexisted=%s, "
             "annihilated=%s, winding_conserved=%s, E_drift=%.4f, motion=%s。%s"
             % (config_label, gamma, steps * DT, SHAPE, n0, nf, coexisted, annihilated,
                winding_conserved, e_drift_frac, motion, notes_extra))
    run_dir = io.write_results(room_id, genesis_yaml, report, integrity, input_vs_output,
                                figures={}, notes=notes)
    print("  wrote %s  [%s g=%.2f] n %d->%d coexist=%s annih=%s Edrift=%.4f motion=%s"
          % (run_dir, config_label, gamma, n0, nf, coexisted, annihilated, e_drift_frac, motion))
    return snaps, phys, report


def main():
    print("=== 渦-依頼A: クリーンな3D渦共生 (媒質依存) ===")
    box = SHAPE[1] * DX
    c = box / 2

    # 逆符号ペア(net-zero): 保存的=共存/並進, 散逸=対消滅。 --- 決定的対照は媒質(gamma)。
    pair = [(c - 2.0, c, +1), (c + 2.0, c, -1)]
    print("\n[逆符号ペア・保存的] gamma=0 → 共存(両生存, 並進)")
    run_and_save("vortex-a-opppair-conservative-seed0001", pair, 0.0, "opposite_pair",
                 notes_extra="保存的媒質: 両渦生存の予言。")
    print("[逆符号ペア・散逸] gamma=0.3 → 対消滅（決定的対照: 媒質が運命を決める）")
    run_and_save("vortex-a-opppair-dissipative-seed0001", pair, 0.30, "opposite_pair",
                 notes_extra="決定的対照: 同一配置でgammaのみ0→0.3。対消滅の予言。")

    # net-zero 4渦(交互符号正方形): 保存的=回転結晶(共生), 散逸=対照。
    d = 4.0
    quad = [(c - d, c - d, +1), (c + d, c - d, -1), (c + d, c + d, +1), (c - d, c + d, -1)]
    print("\n[4渦・保存的] gamma=0 → 回転結晶(4芯すべて生存し回る=共生)")
    run_and_save("vortex-a-crystal-conservative-seed0001", quad, 0.0, "quad_crystal",
                 steps=1600, notes_extra="保存的媒質: 束縛回転結晶(共生)の予言。")
    print("[4渦・散逸] gamma=0.3 → 対照")
    run_and_save("vortex-a-crystal-dissipative-seed0001", quad, 0.30, "quad_crystal",
                 steps=1600, notes_extra="対照: 散逸下の4渦。")

    print("=== 渦-依頼A done ===")


if __name__ == "__main__":
    main()
