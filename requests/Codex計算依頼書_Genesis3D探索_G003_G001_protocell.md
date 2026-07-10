# Codex 計算依頼書 — Genesis 3D 探索: G003 × G001 / protocell 候補

**発行**: Claude（サンドボックス設計・足場）　**宛先**: Codex（正式 3D ソルバの実装・実行）
**リポジトリ**: `Genesis-Room-`（確認用。本体 Aeterna-Genesis とは別、`experiments/e001–e036` は対象外）

---

## §A. 依頼の背景と問い

`docs/ROOM_MODEL.md` §5 は最初の正式 Room 候補として次の2つを挙げている（原文）：

> - **G001：3D 対称性破れ・欠陥形成**（3D TDGL / 複素 Ginzburg–Landau / GPE クエンチ）。期待：Lv1
>   分散/構造因子増大 → Lv2 欠陥/渦線/渦ループ → Lv3 候補 欠陥の移動/再結合。
> - **G003：3D 相分離と流体の共発展**（Model H / Cahn–Hilliard–Navier–Stokes / 反応性
>   phase-field hydrodynamics）。部品接続でなく、同じ自由エネルギー・化学ポテンシャル・応力・
>   流速場から相分離/界面/流れ/変形/輸送が同時に生じる全体場として扱う。最初の統合候補として
>   有望だが G001・G002 より後。

**この依頼の問い**：G001（対称性破れ・位相欠陥形成）と G003（相分離＋流体の共発展）を、それぞれ
**t=0 から中断なく独立に発展させたとき**、以下のいずれかが測定量として自然に現れるかを探索する。

1. G003 単体で、区別可能な領域が周囲から持続し（`docs/EMERGENCE_LEVELS.md` Level 4 の
   「持続する全体性・個体性候補」）、境界・流れ・輸送・変換が同じ場から同時に発生する
   （Level 5「機能の共分化」）——「protocell 候補」と呼べる測定的根拠になるか。
2. G001 の欠陥/渦線形成が、G003 型の相分離場と結合したときに（**Room を分岐させた子 Room として**、
   G001 と G003 を混ぜた新しい genesis condition で。既存の G001/G003 の途中状態へ後付けしない）
   欠陥がドメイン境界の核として働き、Level 4 相当の持続構造を作るか。

**「protocell」は analogy（`docs/PHYSICS_INTEGRITY.md` §9）。** 正式定義（測定量）を満たすまでは
「protocell 候補」「持続構造候補」等の限定名称のみを使うこと。「生命」「細胞」等の強い語は使わない。

---

## §A-2. スコープと非スコープ

- **やること**：G001・G003 を `genesis/` にそれぞれ Room スクリプトとして実装し、正式 3D
  （`full-3d` 相当、`docs/DIMENSION_POLICY.md`）で t=0 から実行し、到達 Level を測定する。
  G001×G003 の結合は**別 genesis condition の子 Room**として、分岐保存（`docs/ROOM_MODEL.md` §3）
  で追加してよい。
- **やらないこと**：完成した膜・境界・分裂位置/時刻を初期条件に置くこと。2D 結果を 3D 結果として
  報告すること。局所 3D／低解像度全体 3D だけで正式 Level 判定を行うこと
  （`docs/DIMENSION_POLICY.md` §1）。
- **2D／局所3D／薄いスラブは事前計算・リスク評価にのみ使ってよい**（`docs/DIMENSION_POLICY.md` §4）。
  正式な到達 Level 判定・§B 報告は本番 3D run から行う。

---

## §A-3. AI が変更してよい範囲（始原側のみ）

`docs/AI_EXPERIMENT_POLICY.md` §2 に従い、次のみ変更可：初期分布・ノイズ強度/相関長・物理定数・
拡散/反応係数・外部流入量/分布・初期対称性・局所相互作用範囲・空間/時間スケール・保存則を満たす
局所法則候補。法則そのものを変える場合は `mutation_type: law_variant` として明示し、より厳しい
監査（物理的由来・対称性・次元整合性・保存則・熱力学整合性・既知極限・第8監査）を通す
（`docs/AI_EXPERIMENT_POLICY.md` §4）。

**変更してはいけないもの**（`docs/ROOM_MODEL.md` §2, `docs/AI_EXPERIMENT_POLICY.md` §3）：完成した
身体形状・膜の完成形・内外指定・分裂位置/時刻・安定形状へ引き寄せる目標関数・求める現象を直接
実行する if 文・中心制御装置・外部 oracle・成功判定コードや監査閾値そのもの。

---

## §B. 返却契約（正典）

各 run について、`common.io.write_results` を使い `results/{run_id}/` に以下を機械的に保存し、
本依頼書への回答として次の4ブロック＋ viz を必ず揃えること。

### B1. Level 測定（`summary.json` の `emergence`）

`common.emergence.compute_level_report` の出力（`docs/EMERGENCE_LEVELS.md` の
`emergence:` schema に沿う）：`reached_level` / `candidate_level` / `detected` / `measured_by`
（画像でなく数値）/ `purity`（`role: E`=純粋創発 か `role: S`=足場付き合成 か。per_object_labels /
external_optimum のいずれかが true なら S）。

Level 4 以上を主張する場合は、`common/emergence.py` の `_DETECTED_KEYS` に合わせて
`persistent_individuality` / `co_differentiation` 等を追記し、根拠となる測定量
（tracked ID lifetime、inside/outside contrast、co-occurring processes の相互依存など）を
`measured_by` に入れること。

### B2. 整合性（`summary.json` の `integrity`、`common.io.integrity_block`）

- 保存量（自由エネルギー/質量など）のドリフト。
- 解像度依存性（最低2段階の格子で比較）。
- 複数 seed（最低3つ）での再現性（到達 Level が一致するか）。
- NaN／無断クリッピングが発生していないこと。

### B3. 入れた vs 出た（`summary.json` の `input_vs_output`、`common.io.input_output_selfcheck`）

第8監査（`docs/PHYSICS_INTEGRITY.md` §6）の4問に機械的に回答：

1. 初期条件に、証明したい量そのもの（またはその関数）が入っていないか？
2. ゲートが、結論の因果を直接 if 判定していないか？
3. ゲートの閾値が、対照（null／線形／ランダム）でも通ってしまわないか？
4. 「創発した」量が、実は別の入力の代数的な言い換えでないか？

**必須の対照（`control_runs` に記録）**：
- 「機構 X（例：G001 の欠陥）が現象 Y（例：ドメイン核形成）を起こす」と主張するなら、X を切った
  対照で Y が消えることを並べる。
- G003 の「持続構造候補」を有界空間の効果でなく本物の閉環として主張するなら、有界空間＋中立
  （選択なし）対照を走らせる（`docs/PHYSICS_INTEGRITY.md` §6「必須の対照」）。

### B4. manifest（`manifest.json`）

`equations` / `solver` / `dt` / `dx` / `grid` / `boundary` / `params` / `seed(s)` / `commit` /
`checksum`（`common.io.checksum_of`）。同じ genesis + seed + mode は同一 checksum になること。

### viz（`figures/`）

可視化と実データの対応を明示する（`docs/PHYSICS_INTEGRITY.md` §2）。可視化データを物理結果として
使わない——figures は補助であり、B1-B4 の数値が正典。

### NOTES.md

出たもの・出なかったもの・EXPLORATION 所見。`docs/PHYSICS_INTEGRITY.md` §2「出なかったものも記録
されている」に対応。7監査（`experiments/*/AUDIT.md` 形式を参考にした自己点検表）を含めてよい。

---

## §C. 昇格段階

`docs/ROOM_MODEL.md` §6 / `docs/AI_EXPERIMENT_POLICY.md` §6 のはしごを飛ばさないこと：

```
IDEA → 2D SCREENED → 2D REPRODUCIBLE → DIMENSION AUDIT PASSED
→ LOCAL 3D PASSED → COARSE GLOBAL 3D PASSED → FULL 3D REPRODUCIBLE
→ PHYSICS AUDIT PASSED → TEMPLATE CANDIDATE → OFFICIAL ROOM
```

このリポジトリの `results/` への保存は「FULL 3D REPRODUCIBLE」以降の段階の記録として扱う。
2D/局所3D の事前計算結果も残してよいが、`mode` フィールドで明示し、正式な到達 Level 判定には
使わない。

---

## §D. 参照

- `docs/AGENTS.md` — 最重要規律・役割分担。
- `docs/EMERGENCE_LEVELS.md` — Level 0–8 と測定指標。
- `docs/PHYSICS_INTEGRITY.md` — 誠実さの規律、role（E/V/S/N/F/Q）、第8監査。
- `docs/ROOM_MODEL.md` — Genesis Room の定義、G001–G003、schema。
- `docs/DIMENSION_POLICY.md` — 2D探索/3D正式/次元移行監査。
- `docs/AI_EXPERIMENT_POLICY.md` — AI の変更可能範囲・昇格段階。
- `genesis/_template.py`, `genesis/README.md` — Room スクリプトの型。
- `common/diagnostics.py`, `common/emergence.py`, `common/io.py` — 測定・保存の共有実装。

不明点や、この依頼書と実装の食い違いがあれば、run を進める前に記録して問い合わせること
（結果を見てから成功条件を変更しない、`docs/PHYSICS_INTEGRITY.md` §3）。
