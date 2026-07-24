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

- `cross_protein_comparison/` —— `dd_idea-run Q8IZL9 P24941 P20794
  --reference Q8IZL9`の出力(詳細は[`../README.jp.md`](../README.jp.md#実例-cdk20-vs-cdk2-vs-mak)参照)。
  CDK20・CDK2・MAKの取得済み配列/AlphaFoldモデル/RCSB構造、蛋白間の配列・
  ポケット保存性マッピング、重ね合わせ済み座標を含む。参照ポケット: 24残基、
  druggability 0.646、位置14・131がCDK2・MAK両方でnon-conservative
  (CDK20選択性を持つ候補残基)。

- `pocket_detection/` —— `dd_idea-search Q8IZL9 -o
  CDK20_inhibition/pocket_detection`の出力(詳細は
  [`../README.jp.md`](../README.jp.md#dd_idea-search-blastベースの入口ツール)参照)。
  CDK20自身はRCSB構造を1件も持たないため、実構造のテンプレートは
  BLASTPで見つけた類似蛋白から取得する——`hits.json`にHomo sapiens限定
  Swiss-Protヒット100件(同一性順)、各行にUniProt family/gene/organism
  メタデータ付き。テンプレート(AlphaFoldモデル+RCSB構造)は
  `dd_idea-search --fetch ACC [ACC ...] -o pocket_detection
  --resolution-cutoff N`でアクセッションごとに選択的に取得する(一括では
  ない)——`hits.json`の各行の`pdb_structures`フィールド参照
  (`null`=未取得)。

## 今後の予定(未着手)

`pocket_detection/hits.json`のどのアクセッションが実際にテンプレート
取得の価値があるか選定する(上位: CDK5 46.0%、CDK3 45.3%、CDK2 43.8%、
CDK7 43.1%、CDK1 43.1%。MAKは35.1%で24位——`dd_idea`のPfam/InterPro
ベース`--discover`では見えないがここでは発見できた)。続いて
`dd_afpocket`のポケット検出/局所拘束MDアンサンブル生成、ドッキング
(`dd_docking`)、QSAR(`dd_chembl`)——実行した際の出力もここに置く
(例: `docking/`、`qsar/`)、各ツール自身のレポジトリには置かない。
