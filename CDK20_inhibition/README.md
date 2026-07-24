[Japanese version](README.jp.md)

# CDK20 inhibitor discovery

Working directory and data for an ongoing CDK20 (human Cyclin-dependent
kinase 20, UniProt `Q8IZL9`) inhibitor-discovery project. Started 2026-07-23
as [`dd_idea`](../)'s own cross-protein-comparison validation run (CDK20 vs.
CDK2 vs. MAK); kept as the real starting point for this project rather than
discarded, since it already established the druggable pocket and how it
differs from the nearest paralogs.

All data for this project lives under this directory -- nowhere else in
this repo, or in any other `dd_*` project's own directory.

## Contents

- `cross_protein_comparison/` -- output of `dd_idea-run Q8IZL9 P24941 P20794
  --reference Q8IZL9` (see [`../README.md`](../README.md#worked-example-cdk20-vs-cdk2-vs-mak)):
  fetched sequences/AlphaFold models/RCSB structures for CDK20, CDK2, and
  MAK, the cross-protein sequence/pocket-conservation mapping, and
  superposed coordinates. Reference pocket: 24 residues, druggability
  0.646, positions 14/131 non-conservative in both CDK2 and MAK (candidate
  CDK20-selectivity residues).

- `pocket_detection/` -- output of `dd_idea-search Q8IZL9 -o
  CDK20_inhibition/pocket_detection` (see [`../README.md`](../README.md#dd_idea-search-blast-based-entry-point)):
  CDK20 has zero RCSB structures of its own, so real structural templates
  for pocket detection/restrained-MD come from BLASTP-similar proteins
  instead -- `hits.json` has 100 Swiss-Prot hits (Homo sapiens, ranked by
  %identity), each with UniProt family/gene/organism metadata. Templates
  (AlphaFold model + RCSB structures) are fetched selectively per
  accession via `dd_idea-search --fetch ACC [ACC ...] -o
  pocket_detection --resolution-cutoff N`, not all at once -- see
  `hits.json`'s `pdb_structures` field per row (`null` = not fetched yet).

## Next steps (not yet done)

Pick which `pocket_detection/hits.json` accessions are actually worth
fetching templates for (top hits: CDK5 46.0%, CDK3 45.3%, CDK2 43.8%,
CDK7 43.1%, CDK1 43.1%; MAK 35.1% at rank 24 -- found here despite being
invisible to `dd_idea`'s Pfam/InterPro-based `--discover`), run
`dd_afpocket`'s pocket detection/restrained-MD ensemble generation,
then docking (`dd_docking`) and QSAR (`dd_chembl`) -- their outputs belong
here too (e.g. `docking/`, `qsar/`), not inside those tools' own repos.
