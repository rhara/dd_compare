[English version](README.md)

# CDK20阻害薬開発

ヒトCDK20(Cyclin-dependent kinase 20、UniProt `Q8IZL9`)の阻害薬開発プロジェクトの
作業ディレクトリ兼データ置き場。2026-07-23に[`dd_idea`](../)自身のクロスプロテイン
比較機能の検証実行(CDK20 vs. CDK2 vs. MAK)として開始したが、既にdruggableな
ポケットと近縁パラログとの違いが判明していたため、破棄せずこのプロジェクトの
実際の出発点として引き継いだ。

このプロジェクトのデータは全てこのディレクトリ配下に置く——このレポジトリの
他の場所や、他の`dd_*`プロジェクト自身のディレクトリには置かない。

## 内容

- `cross_protein_comparison/` —— CDK20・CDK2・MAKの取得済み配列/AlphaFold
  モデル/RCSB構造、蛋白間の配列・ポケット保存性マッピング、重ね合わせ済み
  座標を含む。参照ポケット: 24残基、druggability 0.646、位置14・131が
  CDK2・MAK両方でnon-conservative(CDK20選択性を持つ候補残基)。

- `pocket_detection/` —— CDK20自身はRCSB構造を1件も持たないため、実構造の
  テンプレートはBLASTPで見つけた類似蛋白から取得する。`hits.json`に
  Homo sapiens限定Swiss-Protヒット100件(同一性順)、各行にUniProt
  family/gene/organismメタデータ付き——MAK(同一性35.1%、24位)も含む、
  `dd_idea`のPfam/InterProベース`--discover`では見えない蛋白。各行に
  ChEMBL活性件数(`chembl_targets`、解決済みターゲットごとの
  `n_activities`)も付与済み——CDK20自身はほぼ皆無(活性1件)、詳細は
  下記「分かったこと」参照。構造テンプレート(AlphaFoldモデル+RCSB構造)
  はアクセッションごとに選択的に取得する(一括ではない)——`hits.json`の
  各行の`pdb_structures`フィールド参照(`null`=未取得。現時点では全行
  未取得)。

## 分かったこと

**100件のBLASTヒット間でChEMBLカバレッジに3000倍超の開きがある**
(2026-07-24、`--chembl-activity-all`、pChEMBL値付きbinding assay):
CDK20自身はChEMBL活性データが1件のみ——RCSB構造も0件であることと
整合し、実験的にほぼ未研究のキナーゼであることが分かる。MAK
(ファミリー分類ベースの`--discover`では見えないヒット)は13件。
最もカバレッジが厚いヒットはGSK3B(7448件)、MAPK1/ERK2(6927件)、
MAPK14/p38(6811件)、DYRK1A(6288件)、AURKA(3769件)、CDK2
(3015件、ChEMBL自身の`CHEMBL301`ページと完全一致)。あるヒットが
最大限有用であるためには、CDK20への配列同一性(構造テンプレートの
根拠として)と、十分なSARデータ(`dd_chembl`での後続QSARモデリング用)の
**両方**が必要——`hits.json`の`n_activities`フィールドが後者の
ソート・フィルタ基準になる。

## このディレクトリの再現方法

上記の内容を構築する全コマンドを順番に記載(各コマンドの詳細は
[`../README.jp.md`](../README.jp.md)参照)。どのステップも既にディスクに
キャッシュされているものは再取得せず再利用するので、このリストは中断後の
再開手順としても使える。

2026-07-24、`dd_idea`のツールが当初のアドホックスクリプトから成熟した
ため、クリーンな状態から再構築——`cross_protein_comparison/`は
2026-07-23の初回実行と完全に同一の結果に再構築できた(同じポケット、
同じ24残基、同じdruggability 0.646、同じ蛋白ごとの同一性・保存性の
数値——パイプラインが同一入力に対し決定論的であることを確認)。
`pocket_detection/`のBLASTヒット集合も同一の結果になった(同じ100件、
MAKも同じ24位・35.1%):

```bash
# cross_protein_comparison/ —— クロスプロテイン配列・ポケット比較+構造重ね合わせ
dd_idea-run Q8IZL9 P24941 P20794 -o CDK20_inhibition/cross_protein_comparison --reference Q8IZL9

# pocket_detection/ —— BLASTPベースの類似蛋白テーブル(ダウンロードなし)
dd_idea-search Q8IZL9 -o CDK20_inhibition/pocket_detection

# pocket_detection/ —— 全行についてChEMBLターゲット解決+活性件数カウント
# (軽量——件数のみ、結果は上記「分かったこと」参照)
dd_idea-search --chembl-activity-all -o CDK20_inhibition/pocket_detection
```

**未実行** —— 実際の構造テンプレート取得:

```bash
# pocket_detection/ —— 上記の表を確認した後、特定アクセッションのAlphaFold
# モデル+RCSBテンプレート(2.0Å以下、ツールのデフォルト)を取得。例:
dd_idea-search --fetch P24941 -o CDK20_inhibition/pocket_detection
```

## 今後の予定(未着手)

`pocket_detection/hits.json`の`pct_identity`と`n_activities`の両方を
使って、実際に構造テンプレート取得の価値があるアクセッションを選定する
(候補: CDK2 43.8%/活性3015件、CDK1 43.1%/1488件、AURKA
28.4%/3769件——両軸で高い。MAK 35.1%/13件——クロスプロテイン比較上は
構造的に関連するがChEMBLデータは乏しい)。上記の`--fetch`コマンドを
該当分に実行し、続いて`dd_afpocket`のポケット検出/局所拘束MD
アンサンブル生成、ドッキング(`dd_docking`)、QSAR(`dd_chembl`)——実行
した際の出力もここに置く(例: `docking/`、`qsar/`)、各ツール自身の
レポジトリには置かない。
