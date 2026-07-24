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
  to `dd_idea`'s Pfam/InterPro-based `--discover`. Each row also has its
  ChEMBL bioactivity count (`chembl_targets`, `n_activities` per resolved
  target) -- CDK20 itself has essentially none (1 activity on record);
  see "Findings" below for the well-covered hits. Structural templates
  (AlphaFold model + RCSB structures) are fetched selectively per
  accession, not all at once -- see `hits.json`'s `pdb_structures` field
  per row (`null` = not fetched yet; none fetched so far).

## Findings

**ChEMBL coverage varies by >3000x across the 100 BLAST hits** (2026-07-24,
`--chembl-activity-all`, binding assays with a pChEMBL value): CDK20 itself
has just 1 ChEMBL activity record, consistent with having zero RCSB
structures too -- this is a genuinely under-studied kinase. MAK (the
family-classification-invisible hit `--discover` misses) has 13. The
best-covered hits are GSK3B (7448), MAPK1/ERK2 (6927), MAPK14/p38 (6811),
DYRK1A (6288), AURKA (3769), and CDK2 (3015, matching ChEMBL's own
`CHEMBL301` page exactly). A hit needs *both* decent sequence identity to
CDK20 (for the structural-template rationale) *and* enough SAR data (for
downstream QSAR modeling with `dd_chembl`) to be maximally useful --
`hits.json`'s `n_activities` field is what to sort/filter by for the
latter.

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
# every row (cheap -- counts only, see "Findings" above for the results)
dd_idea-search --chembl-activity-all -o CDK20_inhibition/pocket_detection
```

**Not yet done** -- fetching any actual structural templates:

```bash
# pocket_detection/ -- fetch AlphaFold model + RCSB templates (<=2.0Å, the tool's
# default) for specific accessions once the table above has been reviewed, e.g.:
dd_idea-search --fetch P24941 -o CDK20_inhibition/pocket_detection
```

## Next steps (not yet done)

Using both `pct_identity` and `n_activities` from `pocket_detection/hits.json`,
pick which accessions are actually worth fetching structural templates
for (candidates worth a look: CDK2 43.8%/3015 activities, CDK1 43.1%/1488,
AURKA 28.4%/3769 -- high on both axes; MAK 35.1%/13 -- structurally
relevant per the cross-protein comparison but ChEMBL-poor), run the
`--fetch` command above for those, then `dd_afpocket`'s pocket
detection/restrained-MD ensemble generation, then docking (`dd_docking`)
and QSAR (`dd_chembl`) -- their outputs belong here too (e.g. `docking/`,
`qsar/`), not inside those tools' own repos.
