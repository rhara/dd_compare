"""Cross-protein pairwise sequence alignment: one protein's canonical
sequence (the reference) against another's, indexed by the reference's own
1-based position.

This is the different-proteins counterpart to `dd_seqalign.sequence`'s
`ChainAlignment`/`align_to_canonical`, which aligns every *structure of the
same protein* against that one protein's own canonical sequence. Here there
is no single shared canonical sequence -- each protein has its own -- so a
plain global (not glocal/free-end-gap) alignment is used instead: unlike
same-protein fragments/isoforms, different paralogs can have genuinely
different N-/C-terminal extensions that a free-end-gap alignment would
otherwise silently ignore.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Optional

from Bio.Align import PairwiseAligner, substitution_matrices

_SUBSTITUTION_MATRIX = substitution_matrices.load("BLOSUM62")

_ALIGNER = PairwiseAligner()
_ALIGNER.substitution_matrix = _SUBSTITUTION_MATRIX
_ALIGNER.mode = "global"
_ALIGNER.open_gap_score = -11
_ALIGNER.extend_gap_score = -1


def conservation(a: Optional[str], b: Optional[str]) -> str:
    """'identical' | 'conservative' | 'non-conservative' | 'gap' for a pair
    of one-letter residue codes (either side `None` means a gap -- the
    reference position has no counterpart in the target protein).
    "Conservative" is a positive BLOSUM62 score for the substitution (the
    same scoring matrix the alignment itself uses, rather than a separate
    hand-maintained physicochemical grouping) -- e.g. Phe->Tyr scores +3
    (conservative), Phe->Asp scores -3 (non-conservative)."""
    if a is None or b is None:
        return "gap"
    if a == b:
        return "identical"
    try:
        score = _SUBSTITUTION_MATRIX[a, b]
    except (KeyError, IndexError):
        return "non-conservative"
    return "conservative" if score > 0 else "non-conservative"


def percent_identity(seq_a: str, seq_b: str) -> float:
    """Global BLOSUM62 percent identity between two full sequences, over
    the aligned (non-gap-vs-non-gap) columns only."""
    aln = _ALIGNER.align(seq_a, seq_b)[0]
    s1, s2 = str(aln[0]), str(aln[1])
    aligned = [(x, y) for x, y in zip(s1, s2) if x != "-" and y != "-"]
    if not aligned:
        return 0.0
    matches = sum(1 for x, y in aligned if x == y)
    return 100 * matches / len(aligned)


@dataclass
class ReferenceMapping:
    """A single pairwise global alignment, indexed by the reference
    protein's own 1-based canonical position."""

    reference_accession: str
    target_accession: str
    pct_identity: float
    reference_seq: str
    target_seq: str
    # reference 1-based position -> target 1-based position, or None if the
    # reference position falls opposite a gap (no counterpart in target)
    ref_to_target_pos: Dict[int, Optional[int]] = field(default_factory=dict)

    def target_residue(self, ref_pos: int) -> Optional[str]:
        """The target protein's own residue aligned to `ref_pos` in the
        reference's numbering, or None if there is no counterpart there."""
        pos = self.ref_to_target_pos.get(ref_pos)
        return self.target_seq[pos - 1] if pos else None

    def reference_residue(self, ref_pos: int) -> Optional[str]:
        return self.reference_seq[ref_pos - 1] if 1 <= ref_pos <= len(self.reference_seq) else None


def align_to_reference(
    reference_seq: str, reference_accession: str, target_seq: str, target_accession: str,
) -> ReferenceMapping:
    """Global-align `target_seq` against `reference_seq` and return a
    mapping from every reference position to its counterpart position (or
    None) in the target -- the basis for `pocketmap.py`'s cross-protein
    pocket-residue translation."""
    aln = _ALIGNER.align(reference_seq, target_seq)[0]
    s1, s2 = str(aln[0]), str(aln[1])

    mapping: Dict[int, Optional[int]] = {}
    ref_i = tgt_i = 0
    matches = aligned_cols = 0
    for ca, cb in zip(s1, s2):
        if ca != "-":
            ref_i += 1
        if cb != "-":
            tgt_i += 1
        if ca != "-" and cb != "-":
            mapping[ref_i] = tgt_i
            aligned_cols += 1
            if ca == cb:
                matches += 1
        elif ca != "-":
            mapping[ref_i] = None

    pct = 100 * matches / aligned_cols if aligned_cols else 0.0
    return ReferenceMapping(
        reference_accession=reference_accession, target_accession=target_accession, pct_identity=pct,
        reference_seq=reference_seq, target_seq=target_seq, ref_to_target_pos=mapping,
    )
