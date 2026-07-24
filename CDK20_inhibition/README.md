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

## Next steps (not yet done)

Pocket detection/restrained-MD ensemble generation (`dd_afpocket`), docking
(`dd_docking`), and QSAR (`dd_chembl`) stages for CDK20 have not been run
yet as part of this project -- when they are, their outputs belong here too
(e.g. `pocket_detection/`, `docking/`, `qsar/`), not inside those tools'
own repos.
