"""Generic PDB-chain-to-canonical-UniProt-sequence alignment: given a
downloaded structure, extract each chain's observed sequence, glocal-align
it against a protein's canonical sequence, and pick out which chain is
actually the protein of interest (vs. a bound partner like a cyclin or a
CAK). Vendored (not imported) from `dd_seqalign.sequence`, following this
whole project family's established convention of not depending on any
sibling `dd_*` package at install time.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.Data.IUPACData import protein_letters_3to1
from Bio.PDB import PDBParser
from Bio.PDB.Polypeptide import is_aa

_THREE_TO_ONE = {k.upper(): v for k, v in protein_letters_3to1.items()}

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
