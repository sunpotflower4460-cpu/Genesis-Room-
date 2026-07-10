# Genesis-Room-

**Genesis Room 3D 計算の確認用リポジトリ。** 本体 [Aeterna-Genesis](https://github.com/sunpotflower4460-cpu/Aeterna-Genesis)
とは別リポジトリ（`experiments/e001–e036` は入らない）。

Codex は `requests/` の依頼書に従い、`genesis/_template.py` の型で 3D ソルバを書き、`common/` を
使って `results/` に §B 形式で結果を保存する。規律は `docs/` に置く。

## 構造

```
Genesis-Room-/
├─ docs/          規律ドキュメント（AGENTS / EMERGENCE_LEVELS / PHYSICS_INTEGRITY /
│                  ROOM_MODEL / DIMENSION_POLICY / AI_EXPERIMENT_POLICY）
├─ common/        共有モジュール（diagnostics: Level 測定、emergence: Level レポート、
│                  io: §B 形式の結果書き出し）
├─ genesis/       Genesis Room スクリプト（_template.py が雛形。実際の 3D ソルバは Codex が書く）
├─ requests/      Codex への計算依頼書
├─ results/       結果出力先（run ごとに §B 形式で保存される）
└─ tests/         共有モジュールの最小テスト
```

## 使い方

```bash
pip install -r requirements.txt
python -m pytest tests/                # 共有 diagnostics のテスト
python genesis/_template.py            # 雛形（2D 例）を実行し results/ に §B 出力を作る
```

新しい Room スクリプトの書き方は `genesis/README.md` を参照。**読む順は `docs/AGENTS.md` から。**

## 役割分担

このリポジトリの足場（docs / common / 雛形 / 構造）は Claude Code が作る。**3D の物理ソルバ本体は
Codex が書く**。フル移行（CLI・Observatory アプリ・CI 総入替・room schema の全実装）は本体
Aeterna-Genesis の話であり、ここでは行わない。
