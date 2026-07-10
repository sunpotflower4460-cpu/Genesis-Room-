# genesis/ — 新しい Room スクリプトの書き方

`genesis/_template.py` は雛形（2D 複素 Ginzburg-Landau クエンチの最小例）。**物理ソルバ本体は
Codex が書く**（このリポジトリの足場は雛形の例示のみ）。

## 手順

1. **雛形をコピー**：`genesis/_template.py` を `genesis/<room_name>.py`（例：`genesis/g003_model_h.py`）
   としてコピーする。
2. **genesis condition を書く**：`GENESIS` dict を、実際に計算したい始原条件（fields / 方程式 /
   solver / dimension / grid / dt / dx / boundary / params / initial_state / seed）に差し替える。
   - `docs/ROOM_MODEL.md` §2 の「与えてはならない（完成形）」を確認する。
   - 完成した渦線・欠陥・膜・分裂位置/時刻を初期条件やパラメータに置かない。
   - 3D 正式結果は最初から最後まで 3D で（`docs/DIMENSION_POLICY.md`）。2D/局所3D はあくまで
     事前計算・候補発見であり、正式な到達 Level 判定には使わない。
3. **common を使う**：`common.emergence.compute_level_report` で到達 Level を測定し、Level 4 以上
   を主張する場合は `detected` / `measured_by` に同じキー形式で追記する（`common/emergence.py`
   docstring 参照）。保存量・複数 seed 再現性は `common.io.integrity_block` にまとめる。
4. **§B 保存**：`common.io.write_results` で `results/{run_id}/` に summary.json / manifest.json /
   figures / NOTES.md を書く。`requests/` の依頼書テンプレートの §B が正典。
5. **決定的対照も走らせる**：第8監査（`docs/PHYSICS_INTEGRITY.md` §6）の4問を自己点検し、
   「機構 X を切ると現象 Y が消える」対照や、中立（選択なし）対照を実際に走らせて
   `common.io.input_output_selfcheck` の `control_runs` に記録する。

## 規律（守ること）

- すべての run は時刻 0 から中断なく発展させる。途中でモデルを切り替えない
  （`docs/ROOM_MODEL.md` §1）。
- Level N の後に Level N+1 の機能を"足して"はならない。より深くしたいなら始原条件を変えて
  Level 0 から再実行する（`AGENTS.md`）。
- 過去の Room・結果・失敗を削除しない。失敗（3D で欠陥形成が起きなかった等）も記録として残す
  （`docs/AI_EXPERIMENT_POLICY.md` #9, #12）。
