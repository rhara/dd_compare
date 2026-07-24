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
  ChEMBL活性件数(`chembl_targets`)と解像度2.0Å以下の全RCSB構造
  (`pdb_structures`、`raw_pdb/{gene}/`配下——49gene・930構造、gitignore
  対象・再生成可能でコミットはしていない)も付与済み。これらの組み合わせ
  による優先順位付けは下記「分かったこと」参照。

## 分かったこと

**100件のBLASTヒット間でChEMBLカバレッジに3000倍超の開きがある**
(pChEMBL値付きbinding assay): CDK20自身はChEMBL活性データが1件のみ
——RCSB構造も0件であることと整合し、実験的にほぼ未研究のキナーゼで
あることが分かる。MAK(ファミリー分類ベースの`--discover`では見えない
ヒット)は13件。最もカバレッジが厚いヒットはGSK3B(7448件)、
MAPK1/ERK2(6927件)、MAPK14/p38(6811件)、DYRK1A(6288件)、
AURKA(3769件)、CDK2(3015件、ChEMBL自身の`CHEMBL301`ページと完全
一致)。

**`--rank`(同一性×テンプレート×活性×family、各1〜5の分位数クラスの
積——手法の詳細は[`../README.jp.md`](../README.jp.md#--rank-4つのシグナルを組み合わせる)
参照)では、CDK2とCDK7が満点625で同率トップ**: どちらもCDK20自身と
同じCDC2/CDKX subfamilyに属し、かつ他の全軸でも高スコア(CDK2:
同一性43.8%、テンプレート275件、活性3015件。CDK7: 43.1%、18件、
611件)。次点(スコア400): MAPK14、MAPK1、DYRK1A——ChEMBL/テンプレート
件数は膨大だがCDK20とは別のCMGC subfamily。MAKは本プロジェクトの
最初の実例に登場しfamilyも一致しているにもかかわらず、順位はかなり
下位(68位、スコア40: 同一性クラス4×テンプレートクラス1[2.0Å以下0件]
×活性クラス2[13件]×familyクラス5——ここで足を引っ張っているのは
関連性ではなくデータの乏しさ)——構造的には関連するがデータ不足。
完全なランキング: `pocket_detection/hits_ranked.md`。

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

# pocket_detection/ —— 全行についてChEMBLターゲット解決+活性件数カウント(軽量——件数のみ)
dd_idea-search --chembl-activity-all -o CDK20_inhibition/pocket_detection

# pocket_detection/ —— 全行についてAlphaFoldモデル(シード)+解像度2.0Å以下の
# 全RCSB構造を取得(ツールのデフォルト)——930構造、420MB、数分かかる。
# gitignore対象でコミットはしていない(../README.mdの--fetch-allのコスト
# に関する注記参照)
dd_idea-search --fetch-all -o CDK20_inhibition/pocket_detection

# pocket_detection/ —— 上記4シグナルで全ヒットをランキング(即時・ネット
# ワークアクセスなし)——結果は上記「分かったこと」とhits_ranked.md参照
dd_idea-search --rank -o CDK20_inhibition/pocket_detection --summary-format markdown
```

## 今後の予定(未着手)

`pocket_detection/hits_ranked.md`の上位ヒット(CDK2、CDK7、続いて
MAPK14/MAPK1/DYRK1A)を構造テンプレート/SARデータの優先リストとして、
`dd_afpocket`のポケット検出/局所拘束MDアンサンブル生成、ドッキング
(`dd_docking`)、QSAR(`dd_chembl`)に進む——実行した際の出力もここに
置く(例: `docking/`、`qsar/`)、各ツール自身のレポジトリには置かない。
