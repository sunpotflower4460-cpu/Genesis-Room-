"""snapshots(時系列) -> 到達 Level レポート（docs/EMERGENCE_LEVELS.md の schema に沿う）。"""

import numpy as np

from . import diagnostics as diag

# Keys allowed under emergence.detected -- mirrors the reference emergence.schema.json in the
# main Aeterna-Genesis repo (additionalProperties: false there; kept identical here on purpose).
_DETECTED_KEYS = [
    "difference", "localization", "spontaneous_motion", "circulation",
    "persistent_individuality", "co_differentiation", "self_maintaining_closure",
    "growth_division_inheritance", "selection_open_ended",
]


def compute_level_report(snapshots, kind, per_object_labels=False, external_optimum=False):
    """snapshots(時系列) と kind('cgl'|'model_h'|'protocell'|...) から、EMERGENCE_LEVELS.md の指標を
    計算し、到達 Level と measured 値の dict を返す。role/purity（per_object_labels, external_optimum）
    も含める（純粋 E か足場 S か、EMERGENCE_LEVELS.md「Room への記録」）。

    snapshots: list of dict, 時刻昇順。各要素は少なくとも {"field": array}（2D/3D, 実数 or 複素数）を
      持つ。速度場があれば {"u": [u_x, u_y, (u_z)]} も渡すと Level 3 まで測定する。
    kind: このレポートがどの物理モデルの run かを呼び出し側が記録用に使うラベル（レポート自体には
      含めない -- schema の additionalProperties:false に合わせるため manifest 側で保持する）。
    per_object_labels / external_optimum: docs/EMERGENCE_LEVELS.md Level 5 の「純粋 vs 足場付き」判定。
      いずれか True なら role=S（合成）として記録する。

    このモジュールが機械的に判定するのは Level 0-3（difference/localization/spontaneous_motion/
    circulation）まで。Level 4 以上（persistent_individuality 以降）は Room 固有の追跡（tracked ID、
    load-bearing ablation、分裂検出など）が要るため、呼び出し側が同じ detected/measured_by の形で
    追記すること -- キー名は _DETECTED_KEYS に合わせる。
    """
    if not snapshots:
        raise ValueError("compute_level_report requires at least one snapshot")

    fields = [s["field"] for s in snapshots]
    variances, growth_rate = diag.variance_growth(fields)
    peak_k, prom = diag.structure_factor_peak(fields[-1])
    xi = diag.correlation_length(fields[-1])

    is_complex = np.iscomplexobj(fields[-1])
    defects = diag.winding_defect_count(fields[-1]) if is_complex else 0

    difference = bool(growth_rate > 0 and prom > 1.5 and xi > 0)
    localization = bool(difference and defects > 0)

    detected = {k: False for k in _DETECTED_KEYS}
    detected["difference"] = difference
    detected["localization"] = localization

    measured_by = {
        "variance_growth": round(float(growth_rate), 6),
        "structure_factor_peak_k": round(float(peak_k), 6),
        "structure_factor_prominence": round(float(prom), 6),
        "correlation_length": round(float(xi), 6),
        "defect_count": int(defects),
    }

    reached = 0
    if difference:
        reached = 1
    if localization:
        reached = 2

    if "u" in snapshots[-1]:
        ke = diag.kinetic_energy(snapshots[-1]["u"])
        circ = diag.circulation(snapshots[-1]["u"])
        spontaneous_motion = bool(ke > 0 and circ > 0)
        detected["spontaneous_motion"] = spontaneous_motion
        detected["circulation"] = spontaneous_motion
        measured_by["kinetic_energy"] = round(float(ke), 6)
        measured_by["circulation_proxy"] = round(float(circ), 6)
        if localization and spontaneous_motion:
            reached = 3

    role = "S" if (per_object_labels or external_optimum) else ("E" if reached >= 1 else "F")

    return {
        "reached_level": reached,
        "candidate_level": min(reached + 1, 8),
        "uninterrupted_from_zero": True,
        "level_detected_by_measurement": True,
        "detected": detected,
        "measured_by": measured_by,
        "purity": {
            "per_object_labels": bool(per_object_labels),
            "external_optimum": bool(external_optimum),
            "role": role,
        },
        "natural_emergence": {
            "started_from_time_zero": True,
            "target_shape_seeded": False,
            "runtime_interventions": 0,
            "target_dependent_rules": False,
            "target_dependent_stopping": False,
            "target_dependent_clipping": False,
            "level_detected_by_measurement": True,
        },
    }
