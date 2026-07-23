"""Optional real-PDB-structure overlay: given a protein already fetched as
an AlphaFold model, look up whether it also has real (RCSB) structures and,
if so, pick one and figure out which of its chains/residue numbers
correspond to the canonical UniProt sequence.

This is an *additive* visualization layer only (see README "Why always
AlphaFold, never a PDB structure" and its "real-structure overlay" caveat
below it) -- pocket detection and the cross-protein sequence/pocket mapping
in `pocketmap.py`/`sequence.py` stay anchored on the AlphaFold model
unconditionally, exactly as before. A real PDB entry, when one is picked
here, is shown *alongside* that AlphaFold-based analysis, not instead of it.

RCSB lookup/download (`list_pdb_ids_for_uniprot`, `EntryMetadata`,
`fetch_entry_metadata`, `download_pdb`) and per-chain canonical-sequence
alignment (`ChainSequence`, `extract_chain_sequences`, `ChainAlignment`,
`align_to_canonical`, `pick_target_chain`) are vendored (not imported) from
`dd_seqalign.fetch`/`dd_seqalign.sequence`, following this whole project
family's established convention of not depending on any sibling `dd_*`
package at install time. `dd_seqalign` needs chain-picking to compare every
known structure of *one* protein against itself; here the same logic picks
the right chain out of a real PDB entry that may contain bound partners
(e.g. a cyclin, a CAK) so the *right* protein's numbering is used, and
additionally provides the canonical-position round-trip
(`ChainAlignment.resseq_for_canonical`) needed to place reference pocket-
residue labels onto a PDB entry's own (non-canonical) numbering -- something
`dd_compare`'s AlphaFold-only path never needed, since an AlphaFold model's
numbering already equals canonical UniProt position by construction.
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

from . import pdbio

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"

_THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}


def list_pdb_ids_for_uniprot(accession: str) -> List[str]:
    """Every RCSB PDB entry ID cross-referenced (via SIFTS) to this UniProt
    accession, best-resolution-first (server-side sort) -- a well-studied
    target can have hundreds of entries (e.g. CDK2's 512), and sorting on
    RCSB's side means `select_pdb_structures` only needs to fetch per-entry
    metadata/coordinates for as many of them as it actually scans, instead
    of every one just to rank them."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession.upper(),
            },
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }
    req = urllib.request.Request(
        RCSB_SEARCH, data=json.dumps(query).encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req) as fh:
        body = fh.read()
    if not body:
        return []  # RCSB responds 204 No Content (empty body) when a UniProt accession has zero structures
    return [hit["identifier"] for hit in json.loads(body).get("result_set", [])]


@dataclass
class EntryMetadata:
    pdb_id: str
    method: str
    resolution: Optional[float]
    title: str


def fetch_entry_metadata(pdb_id: str) -> EntryMetadata:
    """Experimental method and resolution (None for methods that don't
    report one, e.g. NMR), straight from RCSB's entry-level summary -- cheap
    (no coordinates), used to rank candidates before downloading any of
    them."""
    with urllib.request.urlopen(RCSB_ENTRY.format(pdb_id=pdb_id.upper())) as fh:
        entry = json.load(fh)
    info = entry.get("rcsb_entry_info", {})
    resolution_list = info.get("resolution_combined") or []
    return EntryMetadata(
        pdb_id=pdb_id.upper(),
        method=info.get("experimental_method", "unknown"),
        resolution=resolution_list[0] if resolution_list else None,
        title=entry.get("struct", {}).get("title", ""),
    )


def download_pdb(pdb_id: str, dest: Path) -> str:
    """Fetch a raw PDB entry from RCSB and return its text contents.
    Cached: skipped if `dest` already exists."""
    dest = Path(dest)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(RCSB_PDB.format(pdb_id=pdb_id.upper()), dest)
    return dest.read_text()


Residue = Tuple[int, str]  # (author resseq, one-letter code)


@dataclass
class ChainSequence:
    chain_id: str
    residues: List[Residue]

    @property
    def sequence(self) -> str:
        return "".join(code for _, code in self.residues)


def extract_chain_sequences(pdb_path: Union[str, Path]) -> Dict[str, ChainSequence]:
    """One `ChainSequence` per protein chain in the structure."""
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure(Path(pdb_path).stem, str(pdb_path))
    model = next(iter(structure))

    chains: Dict[str, ChainSequence] = {}
    for chain in model:
        residues: List[Residue] = []
        for res in chain:
            if not is_aa(res, standard=False):
                continue
            code = _THREE_TO_ONE.get(res.get_resname().strip().upper())
            if code is None:
                continue
            residues.append((res.id[1], code))
        if residues:
            chains[chain.id] = ChainSequence(chain_id=chain.id, residues=residues)
    return chains


@dataclass
class ResidueAlignment:
    canonical_pos: int
    canonical_code: str
    structure_resseq: Optional[int]
    structure_code: Optional[str]
    status: str  # "match" | "mismatch" | "missing"


@dataclass
class ChainAlignment:
    chain_id: str
    residues: List[ResidueAlignment] = field(default_factory=list)
    _reverse: Optional[Dict[int, int]] = field(default=None, repr=False, compare=False)

    @property
    def n_covered(self) -> int:
        return sum(1 for r in self.residues if r.status != "missing")

    @property
    def n_mismatch(self) -> int:
        return sum(1 for r in self.residues if r.status == "mismatch")

    def resseq_for_canonical(self, canonical_pos: int) -> Optional[int]:
        """This chain's own author residue number for a given canonical
        UniProt position, or None if unresolved here -- used to translate a
        reference pocket residue (defined in canonical numbering) into this
        PDB entry's own numbering for label placement."""
        idx = canonical_pos - 1
        if 0 <= idx < len(self.residues):
            return self.residues[idx].structure_resseq
        return None


_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = substitution_matrices.load("BLOSUM62")
_ALIGNER.open_gap_score = -10
_ALIGNER.extend_gap_score = -0.5
_ALIGNER.end_insertion_score = 0.0
_ALIGNER.end_deletion_score = 0.0


def align_to_canonical(chain_seq: ChainSequence, canonical_seq: str) -> ChainAlignment:
    """Glocal-align one chain's observed residues against the canonical
    sequence (free end gaps -- a co-crystal fragment missing its N-/C-
    terminus isn't a real difference from canonical, unlike the *global*
    alignment `sequence.py` uses across different proteins)."""
    alignment = _ALIGNER.align(canonical_seq, chain_seq.sequence)[0]
    result = ChainAlignment(chain_id=chain_seq.chain_id)
    covered = [False] * len(canonical_seq)
    slot: List[Optional[Residue]] = [None] * len(canonical_seq)

    t_blocks, q_blocks = alignment.aligned
    for (t_start, t_end), (q_start, q_end) in zip(t_blocks, q_blocks):
        for i in range(t_end - t_start):
            canon_idx = t_start + i
            chain_idx = q_start + i
            covered[canon_idx] = True
            slot[canon_idx] = chain_seq.residues[chain_idx]

    for canon_idx, canon_code in enumerate(canonical_seq):
        if covered[canon_idx]:
            resseq, code = slot[canon_idx]
            status = "match" if code == canon_code else "mismatch"
            result.residues.append(ResidueAlignment(canon_idx + 1, canon_code, resseq, code, status))
        else:
            result.residues.append(ResidueAlignment(canon_idx + 1, canon_code, None, None, "missing"))
    return result


def pick_target_chain(chain_alignments: Dict[str, ChainAlignment]) -> str:
    """The chain that is actually the protein of interest, ranked by number
    of *matching* residues (not raw coverage) -- picks e.g. the CDK2 chain
    out of a CDK2/cyclin co-crystal, or the CDK2 chain (not a higher-
    coverage-but-mismatching CDK1 chain) out of a CAK assembly."""
    return max(chain_alignments, key=lambda cid: chain_alignments[cid].n_covered - chain_alignments[cid].n_mismatch)


@dataclass
class SelectedPdbStructure:
    accession: str
    pdb_id: str
    resolution: Optional[float]
    chain_id: str
    pdb_path: str
    ligand_resname: Optional[str]  # None if apo
    chain_alignment: ChainAlignment
    n_candidates_scanned: int


def select_pdb_structures(
    accession: str, canonical_seq: str, out_dir: Union[str, Path], *,
    scan_cap: int = 25, min_ligand_atoms: int = 5, max_structures: int = 3,
    resolution_cutoff: float = 2.0, show_progress: bool = True,
) -> List["SelectedPdbStructure"]:
    """Pick up to `max_structures` real PDB structures for `accession`,
    preferring ones with a genuine bound ligand (not water/cryoprotectant/
    cofactor -- see `pdbio.classify_hetero_groups`), one per *distinct*
    ligand (deduplicated by resname, resolution-ranked-first-seen wins) so
    a target crystallized many times with the same fragment doesn't crowd
    out other chemical matter. Falls back to a single best-resolution entry
    overall if none of the (resolution-sorted) first `scan_cap` candidates
    has a ligand at all. Returns `[]` if the accession has no RCSB
    structures. Candidates are scanned cheapest-metadata-first and
    downloaded/parsed for ligand content only as needed, stopping once
    `max_structures` distinct ligands are found or `scan_cap` candidates
    have been scanned -- a well-studied target (e.g. CDK2's 512 entries)
    would otherwise mean downloading hundreds of structures.

    `resolution_cutoff` (Angstrom, lower is better -- 2.0 by default)
    excludes any candidate whose resolution is worse than this or entirely
    unreported (e.g. NMR structures don't have one) *before* downloading
    its coordinates at all, for both the ligand-bound and apo-fallback
    paths -- this is a quality floor, not just a ranking tiebreaker, so a
    low-resolution entry is never picked just because it happened to be
    the first ligand-bound one scanned."""
    out_dir = Path(out_dir)
    pdb_ids = list_pdb_ids_for_uniprot(accession)  # already best-resolution-first (server-side sort)
    if not pdb_ids:
        if show_progress:
            print(f"[pdbstruct] {accession}: no RCSB structures found", flush=True)
        return []

    selections: List[SelectedPdbStructure] = []
    seen_ligands: set = set()
    best_apo: Optional[Tuple[EntryMetadata, str]] = None
    scanned = 0
    for pdb_id in pdb_ids[:scan_cap]:
        if len(selections) >= max_structures:
            break
        scanned += 1
        try:
            meta = fetch_entry_metadata(pdb_id)
        except Exception:
            continue  # a single entry's metadata being unfetchable shouldn't abort the whole scan
        if meta.resolution is None or meta.resolution > resolution_cutoff:
            continue  # doesn't meet the quality bar; skip without downloading coordinates
        dest = out_dir / f"{meta.pdb_id}.pdb"
        try:
            text = download_pdb(meta.pdb_id, dest)
        except Exception:
            continue
        groups = pdbio.classify_hetero_groups(pdbio.collect_hetero_groups(text))
        ligand = pdbio.pick_ligand_of_interest(groups, min_atoms=min_ligand_atoms)
        if ligand is not None:
            if ligand.resname in seen_ligands:
                continue  # already have this ligand from a better-resolution entry; prefer diversity
            seen_ligands.add(ligand.resname)
            selections.append(_finalize_selection(accession, meta, dest, canonical_seq, ligand.resname, scanned))
            if show_progress:
                print(
                    f"[pdbstruct] {accession}: picked {meta.pdb_id} (resolution={meta.resolution}, "
                    f"ligand={ligand.resname}) [{len(selections)}/{max_structures}] after scanning "
                    f"{scanned} candidate(s)", flush=True,
                )
            continue
        if best_apo is None:
            best_apo = (meta, str(dest))

    if selections:
        return selections
    if best_apo is None:
        if show_progress:
            print(
                f"[pdbstruct] {accession}: {scanned} candidate(s) scanned, none met the "
                f"resolution_cutoff={resolution_cutoff}Å bar (or were downloadable)", flush=True,
            )
        return []
    meta, dest = best_apo
    if show_progress:
        print(
            f"[pdbstruct] {accession}: no ligand-bound entry among {scanned} scanned candidate(s), "
            f"falling back to best resolution {meta.pdb_id} ({meta.resolution})", flush=True,
        )
    return [_finalize_selection(accession, meta, Path(dest), canonical_seq, None, scanned)]


def _finalize_selection(
    accession: str, meta: EntryMetadata, pdb_path: Path, canonical_seq: str,
    ligand_resname: Optional[str], scanned: int,
) -> SelectedPdbStructure:
    chains = extract_chain_sequences(pdb_path)
    alignments = {cid: align_to_canonical(cs, canonical_seq) for cid, cs in chains.items()}
    target_chain = pick_target_chain(alignments)
    return SelectedPdbStructure(
        accession=accession, pdb_id=meta.pdb_id, resolution=meta.resolution, chain_id=target_chain,
        pdb_path=str(pdb_path), ligand_resname=ligand_resname, chain_alignment=alignments[target_chain],
        n_candidates_scanned=scanned,
    )


@dataclass
class PdbOverlayResult:
    accession: str
    pdb_id: str
    resolution: Optional[float]
    ligand_resname: Optional[str]
    chain_id: str
    aligned_pdb: str
    rmsd: Optional[float]
    n_aligned_atoms: int
    error: Optional[str] = None


def align_pdb_overlays(
    reference_afdb_path: Union[str, Path], reference_label: str,
    selections: Dict[str, List[SelectedPdbStructure]], out_dir: Union[str, Path], *, show_progress: bool = True,
) -> Dict[str, List[PdbOverlayResult]]:
    """Superpose every protein's selected real PDB structure(s) (`selections`,
    keyed by accession -- proteins with no RCSB structure are simply absent
    or map to an empty list) onto the reference protein's AlphaFold-model
    frame via whole-chain `cealign`, reusing `structalign.align_structures`
    unmodified (it doesn't care whether an input is an AlphaFold model or a
    real structure). The reference's own AFDB file is the fixed anchor; if
    the reference protein itself has selected PDB structure(s) too, those
    get cealigned onto its own AFDB model just like every other protein's
    does -- a real PDB entry's coordinate frame has no relation to
    AlphaFold's, even for the same protein, so each still needs its own
    fit. Each structure gets a unique `structalign` label
    (`f"{acc}_pdb{i}"`) so multiple entries for the same protein don't
    collide, and so `structalign`'s output filenames stay distinct."""
    from . import structalign

    if not any(selections.values()):
        return {acc: [] for acc in selections}

    anchor = structalign.StructureInput(label=reference_label, pdb_path=str(reference_afdb_path), chain_id="A")
    mobiles = [
        structalign.StructureInput(label=f"{acc}_pdb{i}", pdb_path=sel.pdb_path, chain_id=sel.chain_id)
        for acc, sels in selections.items() for i, sel in enumerate(sels)
    ]
    results = structalign.align_structures(
        [anchor] + mobiles, reference_label=reference_label, out_dir=out_dir, show_progress=show_progress,
    )
    by_label = {r.label: r for r in results}

    out: Dict[str, List[PdbOverlayResult]] = {}
    for acc, sels in selections.items():
        entries = []
        for i, sel in enumerate(sels):
            r = by_label.get(f"{acc}_pdb{i}")
            if r is None:
                continue
            entries.append(PdbOverlayResult(
                accession=acc, pdb_id=sel.pdb_id, resolution=sel.resolution, ligand_resname=sel.ligand_resname,
                chain_id=sel.chain_id, aligned_pdb=r.aligned_pdb, rmsd=r.rmsd, n_aligned_atoms=r.n_atoms,
                error=r.error,
            ))
        out[acc] = entries
    return out
