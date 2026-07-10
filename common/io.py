"""依頼書テンプレートの §B 返却契約を機械的に満たす結果書き出し。

results/{run_id}/ に summary.json / manifest.json / figures/*.png / NOTES.md を書く。
§B の4ブロック（Level 測定 / 整合性 / 入れた vs 出た / manifest）+ viz が必ず揃う形にする。
`requests/` の依頼書テンプレートの §B を正典として参照。
"""

import hashlib
import json
import os

import numpy as np


def checksum_of(arrays):
    """SHA-256 checksum over one or more numpy arrays (order-sensitive), hex digest.
    arrays: a single array, or a sequence of arrays (e.g. one per field / seed run).
    """
    if isinstance(arrays, np.ndarray):
        arrays = [arrays]
    h = hashlib.sha256()
    for a in arrays:
        a = np.asarray(a)
        if np.iscomplexobj(a):
            h.update(np.ascontiguousarray(a.real).tobytes())
            h.update(np.ascontiguousarray(a.imag).tobytes())
        else:
            h.update(np.ascontiguousarray(a).tobytes())
    return h.hexdigest()


def integrity_block(conservation_drift, resolutions_result, seed_success, nan_or_clip):
    """docs/PHYSICS_INTEGRITY.md §2 の最低条件のうち機械的に記録できるものをまとめる（§B の
    「整合性」ブロック）。

    conservation_drift: float, 保存量/自由エネルギーの drift（0 に近いほど良い）。
    resolutions_result: dict {resolution_label: measured_value, ...} 解像度依存性の記録。
    seed_success: dict {seed: reached_level, ...} 複数 seed の再現性。
    nan_or_clip: bool, NaN や無断クリッピングが発生したか（False であるべき）。
    """
    seed_levels = list(seed_success.values())
    reproducible = bool(seed_levels) and len(set(seed_levels)) == 1
    return {
        "conservation_drift": float(conservation_drift),
        "resolution_dependence": resolutions_result,
        "seed_reproducibility": seed_success,
        "reproducible_across_seeds": reproducible,
        "nan_or_uncontrolled_clip": bool(nan_or_clip),
        "passed": bool(reproducible and not nan_or_clip),
    }


def input_output_selfcheck(target_encoded_in_initial_condition, gate_encodes_conclusion_causality,
                            gate_passes_null_control, emerged_quantity_is_algebraic_restatement,
                            control_runs):
    """第8監査（docs/PHYSICS_INTEGRITY.md §6）の4問 + 走らせた対照（control_runs）をまとめる
    （§B の「入れた vs 出た」ブロック）。各 boolean 引数は「そう疑われる/そうなってしまった」なら
    True（つまり全部 False が望ましい）。

    control_runs: list of dict, 例 [{"name": "no_X_control", "result": "Y disappears"}]。
    """
    questions = {
        "q1_target_in_initial_condition": bool(target_encoded_in_initial_condition),
        "q2_gate_encodes_conclusion_causality": bool(gate_encodes_conclusion_causality),
        "q3_gate_passes_null_control": bool(gate_passes_null_control),
        "q4_emerged_quantity_is_algebraic_restatement": bool(emerged_quantity_is_algebraic_restatement),
    }
    target_encoded = any(questions.values())
    return {
        "audit8_questions": questions,
        "target_encoded": target_encoded,
        "suggested_role_if_encoded": "Q" if target_encoded else None,
        "control_runs": list(control_runs),
    }


def write_results(run_id, genesis_yaml, emergence_report, integrity, input_vs_output, figures, notes,
                   results_root="results"):
    """§B の4ブロック（Level 測定 / 整合性 / 入れた vs 出た / manifest）+ viz を
    results/{run_id}/ に書き出す。

    genesis_yaml: dict, genesis condition の元
      (equations/solver/dt/dx/grid/boundary/params/seed(s)/commit/checksum)。
    figures: dict {filename: matplotlib Figure}（省略可。空 dict なら figures/ は作らない）。
    notes: str, NOTES.md 本文に足す自由記述（EXPLORATION 所見・出なかったもの）。
    Returns the run directory path.
    """
    run_dir = os.path.join(results_root, run_id)
    os.makedirs(run_dir, exist_ok=True)

    summary = {
        "run_id": run_id,
        "emergence": emergence_report,
        "integrity": integrity,
        "input_vs_output": input_vs_output,
    }
    with open(os.path.join(run_dir, "summary.json"), "w") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)

    manifest = {
        "run_id": run_id,
        "equations": genesis_yaml.get("equations"),
        "solver": genesis_yaml.get("solver"),
        "dt": genesis_yaml.get("dt"),
        "dx": genesis_yaml.get("dx"),
        "grid": genesis_yaml.get("grid"),
        "boundary": genesis_yaml.get("boundary"),
        "params": genesis_yaml.get("params"),
        "seed": genesis_yaml.get("seed"),
        "seeds": genesis_yaml.get("seeds"),
        "commit": genesis_yaml.get("commit"),
        "checksum": genesis_yaml.get("checksum"),
    }
    with open(os.path.join(run_dir, "manifest.json"), "w") as f:
        json.dump(manifest, f, indent=2, ensure_ascii=False)

    if figures:
        fig_dir = os.path.join(run_dir, "figures")
        os.makedirs(fig_dir, exist_ok=True)
        for name, fig in figures.items():
            fig.savefig(os.path.join(fig_dir, name), dpi=150, bbox_inches="tight")

    lines = [
        "# NOTES — %s" % run_id,
        "",
        "## 出たもの（measured_by / detected）",
        "```json",
        json.dumps({"detected": emergence_report.get("detected"),
                    "measured_by": emergence_report.get("measured_by")}, indent=2, ensure_ascii=False),
        "```",
        "",
        "## 出なかったもの・未達",
        "reached_level=%s, candidate_level=%s" % (emergence_report.get("reached_level"),
                                                    emergence_report.get("candidate_level")),
        "",
        "## EXPLORATION 所見",
        notes.strip() if notes else "(なし)",
        "",
    ]
    with open(os.path.join(run_dir, "NOTES.md"), "w") as f:
        f.write("\n".join(lines))

    return run_dir
