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

## 今後の予定(未着手)

ポケット検出/局所拘束MDアンサンブル生成(`dd_afpocket`)、ドッキング
(`dd_docking`)、QSAR(`dd_chembl`)の各段階はこのプロジェクトとしてはまだ
実行していない——実行した際の出力もここに置く(例: `pocket_detection/`、
`docking/`、`qsar/`)、各ツール自身のレポジトリには置かない。
