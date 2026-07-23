"""Detect the reference protein's druggable pocket and map its lining
residues onto every other protein's own numbering.

Because dd_compare's representative structure is always the AlphaFold DB
model (see README "Why always AlphaFold"), every protein's structure residue
number is identical to its own canonical UniProt sequence position: an
AlphaFold model is built directly from the canonical sequence with no gaps,
insertions, or alternate numbering. This means the reference's pocket-lining
residue numbers *are already* reference-canonical positions, and a target
protein's mapped canonical position *is already* that target's own AFDB
model residue number -- no separate structure-numbering round-trip is
needed here, unlike `dd_seqalign.activesite` (which exists specifically to
handle real PDB structures' non-canonical numbering and missing density).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

from .pocket import PocketSelection, find_druggable_pocket
from .sequence import ReferenceMapping, conservation


def detect_reference_pocket(
    reference_pdb: Union[str, Path], work_dir: Union[str, Path], *,
    pocket_rank: int = 1, show_progress: bool = True,
) -> PocketSelection:
    """fpocket-based druggable pocket detection on the reference protein's
    AlphaFold model (no bound ligand needed -- works on any apo structure,
    which every AlphaFold model is)."""
    return find_druggable_pocket(
        Path(reference_pdb), Path(work_dir), pocket_rank=pocket_rank, show_progress=show_progress,
    )


@dataclass
class PocketResidueComparison:
    reference_position: int  # == the reference protein's own AFDB residue number
    reference_residue: str
    target_position: Optional[int]  # == the target protein's own AFDB residue number, if covered
    target_residue: Optional[str]
    conservation: str  # 'identical' | 'conservative' | 'non-conservative' | 'gap'


def compare_pocket(pocket: PocketSelection, mapping: ReferenceMapping) -> List[PocketResidueComparison]:
    """For every reference pocket-lining residue, its counterpart residue/
    position in `mapping`'s target protein (`conservation='gap'` if the
    reference position has no counterpart in the target at all -- e.g. a
    genuine insertion/deletion between the two proteins at that point)."""
    out = []
    for residue in pocket.residues:
        ref_pos = residue.resnum
        ref_code = mapping.reference_residue(ref_pos)
        target_pos = mapping.ref_to_target_pos.get(ref_pos)
        target_code = mapping.target_residue(ref_pos)
        out.append(
            PocketResidueComparison(
                reference_position=ref_pos,
                reference_residue=ref_code or "?",
                target_position=target_pos,
                target_residue=target_code,
                conservation=conservation(ref_code, target_code),
            )
        )
    return out
