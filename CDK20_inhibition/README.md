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

- `cross_protein_comparison/` -- fetched sequences/AlphaFold models/RCSB
  structures for CDK20, CDK2, and MAK, the cross-protein sequence/pocket-
  conservation mapping, and superposed coordinates. Reference pocket: 24
  residues, druggability 0.646, positions 14/131 non-conservative in both
  CDK2 and MAK (candidate CDK20-selectivity residues).

- `pocket_detection/` -- CDK20 has zero RCSB structures of its own, so
  real structural templates for pocket detection/restrained-MD come from
  BLASTP-similar proteins instead. `hits.json` has 100 Swiss-Prot hits
  (Homo sapiens, ranked by %identity), each with UniProt family/gene/
  organism metadata -- including MAK (35.1% identity, rank 24), invisible
  to `dd_idea`'s Pfam/InterPro-based `--discover`. Every row now also has
  its ChEMBL bioactivity count (`chembl_targets`) and every RCSB structure
  <=2.0Å (`pdb_structures`, under `raw_pdb/{gene}/` -- 930 structures
  across 49 genes, gitignored/regenerable, not committed). See "Findings"
  below for how these combine into a priority ranking.

## Findings

**ChEMBL coverage varies by >3000x across the 100 BLAST hits** (binding
assays with a pChEMBL value): CDK20 itself has just 1 ChEMBL activity
record, consistent with having zero RCSB structures too -- this is a
genuinely under-studied kinase. MAK (the family-classification-invisible
hit `--discover` misses) has 13. The best-covered hits are GSK3B (7448),
MAPK1/ERK2 (6927), MAPK14/p38 (6811), DYRK1A (6288), AURKA (3769), and
CDK2 (3015, matching ChEMBL's own `CHEMBL301` page exactly).

**`--rank` (identity x templates x activity x family, each 1-5 quantile
classes multiplied -- see [`../README.md`](../README.md#--rank-combining-the-four-signals)
for the method) puts CDK2 and CDK7 at the top, tied at the maximum score
625**: both are in CDK20's own CDC2/CDKX subfamily *and* score well on
every other axis (CDK2: 43.8% identity, 275 templates, 3015 activities;
CDK7: 43.1%, 18 templates, 611 activities). Next tier (score 400): MAPK14,
MAPK1, DYRK1A -- huge ChEMBL/template counts, but a different CMGC
subfamily than CDK20. MAK, despite being family-matched and part of this
project's original worked example, ranks well down the list (rank 68,
score 40: class 4 identity x class 1 templates [0 at <=2.0Å] x class 2
activity [13] x class 5 family -- exact CDC2/CDKX subfamily match, so the
lowest score component is what actually costs it here, not relatedness)
-- structurally relevant, but data-poor. Full ranked list:
`pocket_detection/hits_ranked.md`.

## Reproducing this directory

Every command to build the contents above, in order, using `dd_idea`'s
console scripts (see [`../README.md`](../README.md) for what each one
does). Re-running any step reuses whatever's already cached on disk
(skips re-fetching), so this list also doubles as how to pick up where a
step left off after an interruption.

2026-07-24, regenerated cleanly from scratch once `dd_idea`'s tooling
matured past its original ad hoc scripts -- `cross_protein_comparison/`
was rebuilt identically to the original 2026-07-23 run (same pocket, same
24 residues, same druggability 0.646, same per-protein identity/
conservation numbers, confirming the pipeline is deterministic given the
same inputs); `pocket_detection/`'s BLAST hit set also came back identical
(same 100 hits, MAK still at rank 24/35.1%):

```bash
# cross_protein_comparison/ -- cross-protein sequence/pocket comparison + structural overlay
dd_idea-run Q8IZL9 P24941 P20794 -o CDK20_inhibition/cross_protein_comparison --reference Q8IZL9

# pocket_detection/ -- BLASTP-based similar-protein table (no downloads yet)
dd_idea-search Q8IZL9 -o CDK20_inhibition/pocket_detection

# pocket_detection/ -- resolve ChEMBL target(s) + count bioactivity records for
# every row (cheap -- counts only)
dd_idea-search --chembl-activity-all -o CDK20_inhibition/pocket_detection

# pocket_detection/ -- fetch AlphaFold model (seed) + every RCSB structure <=2.0Å
# (the tool's default) for every row -- 930 structures, 420MB, several minutes;
# gitignored, not committed (see ../README.md's --fetch-all cost note)
dd_idea-search --fetch-all -o CDK20_inhibition/pocket_detection

# pocket_detection/ -- rank every hit by the four signals above (instant, no
# network access) -- results in "Findings" above and hits_ranked.md
dd_idea-search --rank -o CDK20_inhibition/pocket_detection --summary-format markdown
```

## Next steps (not yet done)

Using `pocket_detection/hits_ranked.md`'s top hits (CDK2, CDK7, then
MAPK14/MAPK1/DYRK1A) as the structural-template/SAR-data priority list,
move to `dd_afpocket`'s pocket detection/restrained-MD ensemble
generation, then docking (`dd_docking`) and QSAR (`dd_chembl`) -- their
outputs belong here too (e.g. `docking/`, `qsar/`), not inside those
tools' own repos.
