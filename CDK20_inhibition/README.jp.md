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
  `dd_idea`のPfam/InterProベース`--discover`では見えない蛋白。
  テンプレート(AlphaFoldモデル+RCSB構造)はアクセッションごとに選択的に
  取得する(一括ではない)——`hits.json`の各行の`pdb_structures`フィールド
  参照(`null`=未取得。現時点では全行未取得)。

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
```

**未実行** —— 実際のテンプレート取得:

```bash
# pocket_detection/ —— 上記の表を確認した後、特定アクセッションのAlphaFold
# モデル+RCSBテンプレート(2.0Å以下、ツールのデフォルト)を取得。例:
dd_idea-search --fetch P24941 -o CDK20_inhibition/pocket_detection
```

## 今後の予定(未着手)

`pocket_detection/hits.json`を確認してどのアクセッションが実際に
テンプレート取得の価値があるか選定する(上位: CDK5 46.0%、CDK3 45.3%、
CDK2 43.8%、CDK7 43.1%、CDK1 43.1%。MAKは35.1%で24位)。上記の
`--fetch`コマンドを該当分に実行し、続いて`dd_afpocket`のポケット検出/
局所拘束MDアンサンブル生成、ドッキング(`dd_docking`)、QSAR
(`dd_chembl`)——実行した際の出力もここに置く(例: `docking/`、`qsar/`)、
各ツール自身のレポジトリには置かない。
