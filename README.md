[Japanese version](README.jp.md)

# dd_compare — Cross-protein sequence alignment meets pocket mapping: two or more UniProt accessions in, active-site divergence hints out

Given two or more UniProt accessions for *different* proteins (paralogs,
off-targets, a kinase family), aligns their canonical sequences, maps one
reference protein's druggable pocket onto every other protein's own
numbering, and superposes their AlphaFold DB models for a visual,
color-coded read on where their active sites actually diverge. Can also
propose candidate similar proteins for a single seed accession (Pfam/
InterPro family membership + sequence-identity ranking), instead of
requiring the accession list up front. Designed as a reusable package, not
tied to any specific target family -- the worked example below (human
CDK20/CDK2/MAK) comes from an actual exploratory drug-discovery session,
but any set of UniProt accessions works.

This is the *different-proteins* counterpart to
[`dd_seqalign`](https://github.com/rhara/dd_seqalign), which compares
every known *structure of one protein* against its own canonical sequence.
Here there is no single shared canonical sequence -- each protein has its
own -- so the reference's pocket, not a fixed set of canonical positions,
is what gets translated across proteins.

- **Fetch (`dd_compare-fetch`)**: canonical sequence + AlphaFold DB model
  for each accession (skips anything already on disk on a re-run).
  `--discover SEED_ACC` instead proposes candidate similar proteins for one
  seed accession (see "Similar-protein discovery" below) and writes
  `candidates.json` -- a proposal only; pick the accessions you actually
  want and pass them to a normal `dd_compare-fetch` call.
- **Align (`dd_compare-align`)**: picks a reference protein (default: the
  first one fetched), runs `fpocket` on its AlphaFold model to detect a
  druggable pocket, pairwise-aligns every other protein's canonical
  sequence against the reference's (Biopython `PairwiseAligner`, BLOSUM62,
  global -- see "Why global, not glocal" below), and maps the reference's
  pocket-lining residues onto each protein's own numbering. Each mapped
  position is classified `identical` / `conservative` / `non-conservative`
  / `gap` (BLOSUM62 score sign, not a hand-maintained substitution-group
  table -- see `dd_compare/sequence.py`). Every protein's AlphaFold model is
  then superposed onto the reference via whole-chain PyMOL `cealign`
  (topology-only, no shared numbering needed -- unlike `dd_seqalign`'s other
  fit mode, this has no cross-protein equivalent of `pair_fit`).
- **Run (`dd_compare-run`)**: fetch + align in one step (explicit accession
  list only -- run `--discover` separately first if you want suggestions).
- **App (`streamlit run app.py -- --report-dir DIR`)**: four tabs --
  Overview (per-protein length/%identity/RMSD table), Active-site
  comparison (the pocket-residue mapping table, colored by conservation),
  Structure overlay (`dd_viewer`'s double-buffered `view3d` component,
  vendored and trimmed -- every protein's AlphaFold model superposed,
  distinctly colored, reference pocket residues highlighted on each),
  Candidates (only shown if `candidates.json` exists in the report
  directory -- the ranked output of a prior `--discover` run).

## Worked example: CDK20 vs. CDK2 vs. MAK

```bash
dd_compare-run Q8IZL9 P24941 P20794 -o data/example_cdk20_cdk2_mak --reference Q8IZL9
streamlit run app.py -- --report-dir data/example_cdk20_cdk2_mak
```

Human CDK20 (Q8IZL9) has no experimental (PDB) structure at all; CDK2
(P24941) has 512; MAK (P20794) has 0 -- exactly the situation that makes
"always compare AlphaFold models, never whichever PDB structure happens to
exist" (see "Why always AlphaFold" below) matter: a fair three-way
comparison isn't possible any other way here. Sequence identity to CDK20:
44.1% (CDK2), 36.5% (MAK). The tool's auto-detected top-ranked pocket (24
residues, druggability 0.65) flags reference positions 14 and 131 as
non-conservative in *both* CDK2 and MAK -- candidate CDK20-distinguishing
pocket residues -- plus position 305 as an outright gap in CDK2 (whose
canonical sequence, at 298 aa, is simply shorter than CDK20's 346).
`data/example_cdk20_cdk2_mak/` is committed as a full worked example
(fetched sequences/structures, `report.json`, superposed coordinates).

## Similar-protein discovery

`--discover` picks the *smallest* Pfam/InterPro family cross-referenced on
the seed's UniProt entry that still has more than one member (same-organism
by default; `--any-organism` lifts that). For CDK20 this resolves to
InterPro `IPR050108` ("CDK" family, 26 human members) rather than the much
broader `Pfam:PF00069` ("Pkinase", 344 members -- "has a kinase domain at
all") or the too-narrow `IPR048002` ("CDK20-like_STKc", 1 member -- just
the seed itself). Every member's canonical sequence is then fetched and
ranked by global BLOSUM62 %identity to the seed.

**Known limitation**: this only finds proteins already classified in the
*same* family as the seed. Running `dd_compare-fetch --discover Q8IZL9`
correctly surfaces CDK2 (44.1% identity) among the top candidates, but
*not* MAK -- MAK is a real, biologically relevant comparison (both
kinases were studied together in this project's own worked example above)
but isn't cross-referenced to the CDK family in InterPro, so it never
enters the candidate pool. Cross-family comparisons like CDK20-vs-MAK need
an explicit accession list.

## Installation

Requires Biopython, pandas, numpy, PyMOL (`pymol2`, importable as a library
-- the conda-forge package is `pymol-open-source`, but its distribution
name as seen by `pip`/`pip show` is `pymol`), and the `fpocket` CLI
(conda-forge only, not on PyPI, invoked as a subprocess). Dedicated conda
env:

```bash
mamba create -n dd_compare -c conda-forge python=3.12 biopython pandas \
    numpy matplotlib py3dmol streamlit pymol-open-source fpocket
conda activate dd_compare

cd dd_compare && pip install --no-deps -e ".[app]"  # [app] adds streamlit/py3Dmol
```

This installs three console commands: `dd_compare-fetch`, `dd_compare-align`,
`dd_compare-run`.

## Usage

```bash
dd_compare-fetch --discover Q8IZL9 -o data/discover   # optional: propose candidates for one seed
dd_compare-run Q8IZL9 P24941 P20794 -o data --reference Q8IZL9
streamlit run app.py -- --report-dir data
```

`--reference` (default: the first accession fetched) picks which protein's
pocket gets detected and which structure everything else is superposed
onto. `--pocket-rank` (default 1, top-ranked by fpocket's Druggability
Score) picks a different pocket on the reference if the top one isn't the
site of interest.

All commands print one line per completed item as it happens; pass
`--no-progress` to suppress this and only print the final summary.

## Design notes

- **Vendored, not imported, from sibling `dd_*` projects**: `pdbio.py`
  (HETATM classification), `pocket.py` (the `fpocket` wrapper), and the
  `cealign` half of `structalign.py` are copied in (from `dd_seqalign`,
  itself originally from `dd_prep`/`dd_afpocket`) rather than depending on
  any sibling package at install time -- the whole `dd_*` family's
  established convention, so each project stays installable on its own
  even though several envs' package pins have already diverged in practice.
  `dd_seqalign`'s other structural-fit mode, `cmd.pair_fit` on a known 1:1
  residue correspondence, is *not* carried over: it relies on both sides
  already sharing canonical UniProt positions, which only makes sense for
  multiple structures of the *same* protein.
- **Why always AlphaFold, never a PDB structure**: picking "the best
  available PDB structure, falling back to AlphaFold" would make quality
  and ligand-bound state vary protein-by-protein (see the CDK20/CDK2/MAK
  numbers above -- 0 vs. 512 vs. 0 real structures), confounding "is this
  active-site difference real biology, or an artifact of comparing a real
  co-crystal against an apo model?" Using the AlphaFold model uniformly
  keeps every comparison apples-to-apples. A useful side effect: an
  AlphaFold model's residue numbering is always identical to its own
  canonical UniProt sequence position (no gaps, insertions, or alternate
  numbering), so a pocket residue detected on the reference is already
  expressed in the same coordinate system the cross-protein sequence
  alignment itself uses -- no separate structure-numbering round-trip is
  needed here, unlike `dd_seqalign.activesite` (which exists specifically
  to handle real PDB structures' non-canonical numbering).
- **Why global, not glocal, sequence alignment**: `dd_seqalign` aligns
  fragments/isoforms of the *same* protein against its own canonical
  sequence with free end gaps (glocal), since a co-crystal fragment missing
  its N-/C-terminus isn't a real biological difference. Different paralogs
  genuinely can have different N-/C-terminal extensions (compare CDK20's
  346 aa to MAK's 623), so a glocal alignment would silently discount real
  sequence differences at the ends -- a plain global alignment is used
  instead.
- **Conservation via BLOSUM62 score sign, not a hand-maintained group
  table**: a substitution is "conservative" if its BLOSUM62 score is
  positive (e.g. Phe->Tyr: +3) and "non-conservative" otherwise (e.g.
  Phe->Asp: -3) -- reusing the same scoring matrix the alignment itself is
  built on, rather than maintaining a separate physicochemical grouping
  that could disagree with it.
- **No rdkit, no C++, no GPU**: this project has no small-molecule handling
  (unlike `dd_seqalign`, which vendors `dd_viewer` whole and therefore
  carries `rdkit` along for a single-receptor scene builder dd_compare
  never uses -- only the plain string-patching + component pieces of
  `dd_viewer` are vendored here, see `dd_compare/viewer3d/__init__.py`).
  Its own work is REST fetches and sequence alignment (both cheap) plus
  delegating anything heavier to already-compiled tools (`fpocket`, PyMOL's
  C++ core) -- there's no hot numeric loop of this project's own worth
  porting to C++, and no GPU-parallelizable workload like `dd_afpocket`'s MD
  sampling or `dd_docking`'s docking runs.

## Known limitations

- **Pfam/InterPro-based discovery only finds same-family proteins** (see
  "Similar-protein discovery" above) -- a functionally related but
  differently classified protein (MAK, relative to CDK20/CDK2) needs an
  explicit accession list.
- **fpocket's auto-detected pocket can miss residues that only contribute
  via backbone atoms**, not side chains -- e.g. a kinase hinge residue
  whose side chain points away from the ATP pocket while its backbone N-H/
  C=O still hydrogen-bonds the adenine. Such a residue can be a genuine,
  literature-known selectivity determinant without ever appearing in
  `--pocket-rank 1`'s lining-residue list. `--pocket-rank` can be used to
  try a different fpocket-ranked pocket, but there is currently no way to
  force a specific residue into the comparison if fpocket's geometric
  pocket detection doesn't select it.
- **Always the AlphaFold model, never a real ligand-bound conformation**
  (a deliberate choice, see "Design notes" above) -- induced-fit pocket
  changes upon ligand binding are out of scope for this comparison.
