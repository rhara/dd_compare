[Japanese version](README.jp.md)

# dd_idea — Cross-protein sequence alignment meets pocket mapping: two or more UniProt accessions in, active-site divergence hints out

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

- **Fetch (`dd_idea-fetch`)**: canonical sequence + AlphaFold DB model
  for each accession (skips anything already on disk on a re-run), plus
  (unless `--no-pdb-overlay`) looks up and selects each protein's real RCSB
  structures too (see "Real-structure overlay" below) -- all genuinely
  fetch-time network work, all cached, recorded in `manifest.json`.
  `--discover SEED_ACC` instead proposes candidate similar proteins for one
  seed accession (see "Similar-protein discovery" below) and writes
  `candidates.json` -- a proposal only; pick the accessions you actually
  want and pass them to a normal `dd_idea-fetch` call.
- **Align (`dd_idea-align`)**: picks a reference protein (default: the
  first one fetched), runs `fpocket` on its AlphaFold model to detect a
  druggable pocket, pairwise-aligns every other protein's canonical
  sequence against the reference's (Biopython `PairwiseAligner`, BLOSUM62,
  global -- see "Why global, not glocal" below), and maps the reference's
  pocket-lining residues onto each protein's own numbering. Each mapped
  position is classified `identical` / `conservative` / `non-conservative`
  / `gap` (BLOSUM62 score sign, not a hand-maintained substitution-group
  table -- see `dd_idea/sequence.py`). Every protein's AlphaFold model,
  and every real structure `dd_idea-fetch` already selected for it, is
  then superposed onto the reference via whole-chain PyMOL `cealign`
  (topology-only, no shared numbering needed -- unlike `dd_seqalign`'s other
  fit mode, this has no cross-protein equivalent of `pair_fit`). This step
  makes no RCSB network calls itself -- re-running align (e.g. to try a
  different `--reference`/`--pocket-rank`) only reuses what fetch already
  cached.
- **Run (`dd_idea-run`)**: fetch + align in one step (explicit accession
  list only -- run `--discover` separately first if you want suggestions).
- **App (`streamlit run app.py -- --report-dir DIR`)**: four tabs --
  Overview (per-protein length/%identity/RMSD table), Active-site
  comparison (the pocket-residue mapping table, colored by conservation),
  Structure overlay (a double-buffered `view3d` component vendored and
  trimmed from the now-retired `dd_viewer` project -- every protein's
  AlphaFold model superposed,
  distinctly colored, reference pocket residues highlighted on each, with
  an optional text label per pocket residue -- toggle both the highlight
  and the labels from the sidebar, independent of which proteins are
  shown), Candidates (only shown if `candidates.json` exists in the report
  directory -- the ranked output of a prior `--discover` run). Every real
  RCSB structure a protein has (see "Real-structure overlay" below) gets
  its own checkbox, nested under that protein's own in the "Proteins to
  show" sidebar list -- shown/hidden independently of its AlphaFold model
  and of each other, drawn as a solid, distinctly-colored cartoon with its
  own bound ligand as sticks in the same color (see the color swatches in
  the caption above the 3D view).

## Worked example: CDK20 vs. CDK2 vs. MAK

```bash
dd_idea-run Q8IZL9 P24941 P20794 -o CDK20_inhibition/cross_protein_comparison --reference Q8IZL9
streamlit run app.py -- --report-dir CDK20_inhibition/cross_protein_comparison
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
`CDK20_inhibition/cross_protein_comparison/` is committed as a full worked
example (fetched sequences/structures, `report.json`, superposed
coordinates) -- part of an ongoing CDK20 inhibitor-discovery project that
this tool's own validation run kicked off; see `CDK20_inhibition/README.md`.
CDK2 also has three real, ligand-bound RCSB structures (`6Q4G`/`HJK` 0.98 Å,
`6Q49`/`HGQ` 1.00 Å, `6Q4H`/`HGH` 1.00 Å) picked automatically by the
real-structure overlay (see below) -- the app draws each as its own
solid, distinctly-colored cartoon alongside CDK2's AlphaFold model, so
the reference's mapped pocket residues can be checked against where
three different real, known ligands actually sit.

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
*same* family as the seed. Running `dd_idea-fetch --discover Q8IZL9`
correctly surfaces CDK2 (44.1% identity) among the top candidates, but
*not* MAK -- MAK is a real, biologically relevant comparison (both
kinases were studied together in this project's own worked example above)
but isn't cross-referenced to the CDK family in InterPro, so it never
enters the candidate pool. Cross-family comparisons like CDK20-vs-MAK need
an explicit accession list -- or see `dd_idea-search` below, which finds
MAK directly since it doesn't depend on family classification at all.

## `dd_idea-search`: BLAST-based entry point

A different, complementary way to find similar proteins and pull together
pocket-detection inputs for a new target -- BLASTP against Swiss-Prot
instead of Pfam/InterPro family membership, so it isn't subject to
`--discover`'s known limitation above (it finds MAK for CDK20, at 35.1%
local identity, rank 24 among Homo sapiens Swiss-Prot hits). Takes a
UniProt accession, a ChEMBL target ID (e.g. `CHEMBL301`), or a raw
amino-acid sequence pasted directly -- useful when starting from a ChEMBL
target of interest or a sequence that isn't in UniProt yet.

Three explicit steps, deliberately not one: build the table first (fast,
no downloads), review it, then fetch AlphaFold models/RCSB structures
and/or count ChEMBL bioactivity data only for the accessions actually
worth it.

```bash
dd_idea-search Q8IZL9 -o CDK20_inhibition/pocket_detection
# -> prints a full table (Family/Gene/Organism/Length/%Id/E-value), no downloads

dd_idea-search --fetch P24941 P20794 -o CDK20_inhibition/pocket_detection --resolution-cutoff 2.5
# -> AlphaFold model (if it's the seed) + every RCSB structure <= 2.5Å for just those two

dd_idea-search --fetch-all -o CDK20_inhibition/pocket_detection --resolution-cutoff 2.5
# -> same, for every row in the table (can mean hundreds of structures for a well-studied hit -- see below)

dd_idea-search --chembl-activity-all -o CDK20_inhibition/pocket_detection
# -> resolves each row's ChEMBL SINGLE PROTEIN target(s) and counts binding-assay
#    activities with a pChEMBL value (same filter dd_chembl itself uses for QSAR
#    training data) -- cheap (counts only, one request per target, no bioactivity
#    data actually downloaded), so --all is reasonable here unlike --fetch-all
```

A well-studied hit can have hundreds of RCSB entries (CDK2 alone: 512
total, ~450 at <=2.5Å) -- `--fetch`ing everything indiscriminately isn't
usually what you want; `--fetch` a short, deliberate list instead of
reaching for `--fetch-all`. `--chembl-activity` has no equivalent cost
concern (a handful of small REST calls per accession), so `--chembl-activity-all`
is fine as a default first pass -- e.g. for CDK20's 100 BLAST hits, this
surfaced a >3000x spread in ChEMBL coverage (GSK3B: 7448 activities;
CDK20 itself: 1) worth knowing about before deciding which hits are
useful pocket-detection templates *and* have enough SAR data for follow-up
QSAR work.

`--fetch`ed RCSB structures land under `raw_pdb/{gene}/` (e.g.
`raw_pdb/CDK2/6Q4G.pdb`), not flat in `raw_pdb/` -- a multi-accession
`--fetch-all` can mean thousands of files from dozens of unrelated
proteins otherwise. Falls back to the accession itself when a row has no
gene name on record (e.g. a raw-sequence seed with no UniProt entry). A
PDB entry cross-referenced to more than one accession's UniProt record
(e.g. a multi-chain structure) gets its own copy in each owning gene's
subdirectory, so every gene folder stays self-contained.

### `--rank`: combining the four signals

Once a hit set has identity, template, activity, and family data,
`--rank` is a fourth mode -- pure computation over the already-gathered
`hits.json`, no network access -- that combines all four into a single
ranking:

```bash
dd_idea-search --rank -o CDK20_inhibition/pocket_detection
```

The four signals live on very different scales (ChEMBL activity counts
span 0-7000+; RCSB template counts are zero for roughly half of any
typical hit set; %identity is a compact, roughly continuous 20-50% band),
so `dd_idea.rank` converts each to an ordinal 1-5 class *before*
combining, rather than multiplying raw magnitudes (which would let
whichever metric happens to have the largest numbers dominate regardless
of what it means biologically):

- **Identity class** -- quantile bin of `%identity` among the hit set.
- **Template/activity class** -- zero always gets class 1 on its own
  (about half of CDK20's hits have zero RCSB structures, so folding that
  into the bottom quantile would compress every real nonzero count into
  one or two classes); nonzero values are then quantile-binned among just
  the other nonzero values.
- **Family class** -- how deep a hit's UniProt family hierarchy
  (superfamily -> family -> subfamily) matches the seed's, from the top:
  exact subfamily match -> top class; same family, different subfamily
  -> one below; only the superfamily matches -> the middle class.

The composite score is the product of the four classes (1-625). For
CDK20, this correctly separates "close paralog with real data" (CDK2 and
CDK7, both scoring the maximum 625: identical CDC2/CDKX subfamily, strong
identity, dozens of templates, hundreds+ of activities) from "well-studied
kinase but a more distant relative" (MAPK14/MAPK1/DYRK1A, score 400 --
huge ChEMBL/template counts, but a different CMGC subfamily) and from
"structurally interesting but data-poor" (e.g. MAK, family-matched and
part of the original worked example, but with only 13 ChEMBL activities
and no RCSB structures of its own).

### Recommended workflow: rank before you fetch

`--rank` doesn't require `--fetch`/`--fetch-all` to have run first --
a row whose `pdb_structures` is still `None` (not yet checked) scores
template class 1, the same as a row that *was* checked and had zero, so
it never breaks the ranking, it just can't yet distinguish on that one
axis. Since identity comes free with the initial search and ChEMBL
activity counts are cheap (a couple of small REST calls per accession),
ranking on those two plus family -- *before* spending any time on RCSB --
is enough to see who's actually worth fetching structures for:

```bash
dd_idea-search Q8IZL9 -o CDK20_inhibition/pocket_detection
dd_idea-search --chembl-activity-all -o CDK20_inhibition/pocket_detection
dd_idea-search --rank -o CDK20_inhibition/pocket_detection --summary-format markdown
# -> review hits_ranked.md, no PDB downloads have happened yet

dd_idea-search --fetch P24941 P50613 Q16539 P28482 Q13627 Q00535 P06493 Q00534 Q13164 P49841 \
    P53779 P45984 Q9Y463 P45983 Q14004 P68400 P49759 O43781 Q15759 P27361 \
    -o CDK20_inhibition/pocket_detection
# -> the top 20 accessions read straight off hits_ranked.md, pasted in by hand.
#    Fetches structures for only those, instead of --fetch-all's every hit
#    (CDK20's 99 hits: --fetch-all took several minutes and 420MB just for
#    RCSB metadata scans/downloads; a top-20 shortlist is a fraction of that)

dd_idea-search --rank -o CDK20_inhibition/pocket_detection
# -> re-run to fold the now-known template counts into the ranking too
```

`--fetch-all` is still there for when you actually want everything (e.g.
a final, exhaustive pass once a target's shortlist is settled) -- it's
just not the first thing to reach for.

## Installation

Requires Biopython, pandas, numpy, PyMOL (`pymol2`, importable as a library
-- the conda-forge package is `pymol-open-source`, but its distribution
name as seen by `pip`/`pip show` is `pymol`), and the `fpocket` CLI
(conda-forge only, not on PyPI, invoked as a subprocess). Dedicated conda
env:

```bash
mamba create -n dd_idea -c conda-forge python=3.12 biopython pandas \
    numpy matplotlib py3dmol streamlit pymol-open-source fpocket
conda activate dd_idea

cd dd_idea && pip install --no-deps -e ".[app]"  # [app] adds streamlit/py3Dmol
```

This installs three console commands: `dd_idea-fetch`, `dd_idea-align`,
`dd_idea-run`.

## Usage

```bash
dd_idea-fetch --discover Q8IZL9 -o data/discover   # optional: propose candidates for one seed
dd_idea-run Q8IZL9 P24941 P20794 -o data --reference Q8IZL9
streamlit run app.py -- --report-dir data
```

`--reference` (default: the first accession fetched) picks which protein's
pocket gets detected and which structure everything else is superposed
onto. `--pocket-rank` (default 1, top-ranked by fpocket's Druggability
Score) picks a different pocket on the reference if the top one isn't the
site of interest.

`dd_idea-fetch`/`-run` additionally take (see "Real-structure overlay"
below for what they control): `--no-pdb-overlay` skips the real-RCSB-
structure lookup entirely (fetch/align only against AlphaFold models, as
in earlier versions); `--pdb-max-structures N` (default 3) caps how many
distinct ligand-bound real structures are kept per protein;
`--pdb-resolution-cutoff N` (default 2.0, Angstrom) excludes any candidate
worse than this resolution -- or with no reported resolution at all, e.g.
NMR structures -- before it's even downloaded, for both the ligand-bound
and best-resolution-fallback paths; `--pdb-scan-cap N` (default 25) caps
how many resolution-ranked candidates get checked for a bound ligand, at
most, before giving up on finding `--pdb-max-structures` of them and
falling back to the single best-resolution one, for a target with
hundreds of structures. `dd_idea-align` doesn't take any of these --
it never makes RCSB network calls itself, only reusing whatever
`dd_idea-fetch` already cached (see below).

All commands print one line per completed item as it happens; pass
`--no-progress` to suppress this and only print the final summary.

## Real-structure overlay

`dd_idea-fetch`/`-run` look up, for every protein, whether it has any
real RCSB structures at or better than `--pdb-resolution-cutoff` and -- if
so -- select up to `--pdb-max-structures` of them, one per *distinct*
bound ligand (not water/cryoprotectant/cofactor; the best-resolution entry
for each ligand wins if the same ligand shows up in more than one
candidate), falling back to a single best-resolution entry (still subject
to the same cutoff) if none of the scanned candidates has a ligand at all
-- recorded in `manifest.json`. This selection is genuinely fetch-time
network work, cached exactly like the canonical-sequence/AlphaFold-model
downloads: **want more real structures for a protein you already
fetched?** just re-run `dd_idea-fetch` for the same `-o` directory with
a larger `--pdb-max-structures` (and/or looser `--pdb-resolution-cutoff`)
-- already-downloaded entries are reused, only the additional candidates
needed get fetched:

```bash
dd_idea-fetch Q8IZL9 P24941 P20794 -o data --pdb-max-structures 6
```

`dd_idea-align` then superposes every already-selected real structure
onto the reference alongside that protein's AlphaFold model, with no RCSB
calls of its own. This is purely an *additional* visualization layer:
pocket detection and the cross-protein sequence/pocket mapping stay
anchored on the AlphaFold model unconditionally, exactly as before -- see
"Why always AlphaFold" below, which still holds for the actual comparison.
`report.json` records the picks per protein as a list,
`pdb_structures: [{pdb_id, resolution, ligand_resname, ...}, ...]`; a
protein with no RCSB structures at all (e.g. CDK20 itself) simply has
`pdb_structures: []`.

**Viewing it**: the committed CDK20/CDK2/MAK example already includes a
real-structure overlay for CDK2 -- three distinct ligand-bound structures
-- no extra fetch needed, just open the app against it:

```bash
streamlit run app.py -- --report-dir CDK20_inhibition/cross_protein_comparison
```

Open the **Structure overlay** tab. In the sidebar's "Proteins to show"
list, each real structure gets its own indented checkbox under its
protein's (e.g. `6Q4G (HJK, 0.98Å)` indented under `P24941`), independent
of that protein's own AlphaFold checkbox and of its other real structures
-- check as many or as few as you want to compare at once. Each shown
structure is drawn as a solid cartoon in its own distinct color (cartoon
and bound-ligand sticks alike, cycled from a palette with no black/gray
-- achromatic colors are hard to read), a colored swatch in the caption
above the 3D view identifies which color is which PDB ID/ligand -- a
direct visual check of whether the reference's
mapped pocket residues line up with where several different real, known
ligands actually sit. The caption also names which proteins (CDK20, MAK)
had no RCSB structure at all. The 3D view's camera position is preserved
across every checkbox toggle (it only resets on "Reset view"), so
switching which structures are shown never disrupts a comparison you're
mid-way through setting up.

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
  multiple structures of the *same* protein. `pdbstruct.py` (the
  real-structure-overlay lookup, see above) additionally vendors
  `dd_seqalign.fetch`'s RCSB search/download functions and
  `dd_seqalign.sequence`'s per-chain canonical-sequence alignment
  (`pick_target_chain`, `align_to_canonical`) -- needed here to pick the
  right chain out of a real PDB entry that may include a bound partner
  (e.g. a cyclin) and to translate a reference pocket residue (defined in
  canonical numbering) into that entry's own, non-canonical residue
  numbers for label placement, something the AlphaFold-only path never
  needed (see below).
- **Why always AlphaFold, never a PDB structure -- for the comparison
  itself**: picking "the best available PDB structure, falling back to
  AlphaFold" would make quality and ligand-bound state vary
  protein-by-protein (see the CDK20/CDK2/MAK numbers above -- 0 vs. 512
  vs. 0 real structures), confounding "is this active-site difference real
  biology, or an artifact of comparing a real co-crystal against an apo
  model?" Using the AlphaFold model uniformly keeps pocket detection and
  the cross-protein sequence/pocket mapping apples-to-apples -- this still
  holds unconditionally even with the real-structure overlay above, which
  is display-only and never feeds back into either. A useful side effect:
  an AlphaFold model's residue numbering is always identical to its own
  canonical UniProt sequence position (no gaps, insertions, or alternate
  numbering), so a pocket residue detected on the reference is already
  expressed in the same coordinate system the cross-protein sequence
  alignment itself uses -- no separate structure-numbering round-trip is
  needed for that part, unlike `dd_seqalign.activesite` (which exists
  specifically to handle real PDB structures' non-canonical numbering) and
  unlike `pdbstruct.py`'s own overlay-label placement, which does need
  exactly that round-trip since it's translating onto a real structure's
  numbering.
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
  (unlike `dd_seqalign`, which vendors the now-retired `dd_viewer` project
  whole and therefore carries `rdkit` along for a single-receptor scene
  builder dd_idea never uses -- only the plain string-patching +
  component pieces of `dd_viewer` are vendored here, directly from that
  now-retired project, see `dd_idea/viewer3d/__init__.py`). Its own
  work is REST fetches and sequence alignment (both cheap) plus
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
  for the actual comparison (a deliberate choice, see "Design notes"
  above) -- induced-fit pocket changes upon ligand binding are out of
  scope there; the real-structure overlay is display-only.
- **Real-structure overlay picks up to `--pdb-max-structures` entries per
  protein** (default 3, one per distinct bound ligand) among the
  best-resolution `--pdb-scan-cap` candidates at or better than
  `--pdb-resolution-cutoff` -- it does not scan every structure for a
  target with hundreds of them (e.g. CDK2's 512), so a ligand-bound entry
  outside that scanned window could be missed, and it does not show every
  ligand a target has ever been crystallized with, just the ones picked.
  Deduplication is by exact ligand resname, so two genuinely different
  ligands that happen to share the same three-letter RCSB code (rare, but
  not impossible) would be treated as one. A target with no structures at
  or better than the resolution cutoff gets no real-structure overlay at
  all, the same as one with no RCSB structures whatsoever.
