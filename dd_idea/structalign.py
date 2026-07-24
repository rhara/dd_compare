"""Structural superposition of every protein's AlphaFold model onto one
reference, in one shared PyMOL session (`pymol2`, not the classic global
`pymol`/`cmd` singleton -- keeps this importable as a library without
launching a GUI or fighting other code over global state).

Vendored (cealign path only) from `dd_seqalign.structalign`: PyMOL's CE
algorithm (`cmd.cealign`) finds the best structural superposition on its own
topology, with no shared residue numbering assumed -- exactly what's needed
here, since every input is a *different* protein (not multiple structures of
the same one). `dd_seqalign`'s other path, `cmd.pair_fit` on a known 1:1
residue correspondence, relies on both sides already sharing canonical
UniProt positions, which only makes sense for structures of the *same*
protein -- it has no cross-protein equivalent and is not carried over here.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Union


@dataclass
class StructureInput:
    label: str  # the UniProt accession
    pdb_path: str
    chain_id: str


@dataclass
class AlignmentResult:
    label: str
    reference_label: str
    rmsd: Optional[float]
    n_atoms: int  # cealign's own aligned-residue count
    aligned_pdb: str
    error: Optional[str] = None  # set instead of raising when this one structure can't be fit


def align_structures(
    structures: Sequence[StructureInput], reference_label: str, out_dir: Union[str, Path], *,
    show_progress: bool = True,
) -> List[AlignmentResult]:
    """Superpose every structure in `structures` onto `reference_label` via
    whole-chain `cealign` and save each (including the unmoved reference,
    for a consistent output set) to `out_dir/{label}_aligned.pdb`.

    `show_progress` prints one line per completed fit (`print(...,
    flush=True)`) -- PyMOL's own `cmd.cealign` return value already reports
    RMSD/alignment length, so no separate progress source is needed beyond
    labeling which structure it was for.
    """
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    by_label = {s.label: s for s in structures}
    if reference_label not in by_label:
        raise ValueError(f"reference_label {reference_label!r} not found in structures")
    reference = by_label[reference_label]

    import pymol2

    results: List[AlignmentResult] = []
    with pymol2.PyMOL() as session:
        cmd = session.cmd
        for s in structures:
            cmd.load(s.pdb_path, s.label)

        ref_out = out_dir / f"{reference.label}_aligned.pdb"
        cmd.save(str(ref_out), reference.label)
        results.append(AlignmentResult(reference.label, reference.label, 0.0, 0, str(ref_out)))
        if show_progress:
            print(f"[structalign] {reference.label}: reference (unmoved) -> {ref_out.name}", flush=True)

        mobiles = [s for s in structures if s.label != reference_label]
        for i, s in enumerate(mobiles, start=1):
            try:
                result = cmd.cealign(
                    f"{reference.label} and chain {reference.chain_id} and alt ''+A",
                    f"{s.label} and chain {s.chain_id} and alt ''+A",
                )
            except Exception as e:
                # One structure genuinely not comparable (e.g. cealign finds
                # no usable structural correspondence at all) shouldn't
                # abort the fit for every other structure in the batch.
                results.append(AlignmentResult(s.label, reference_label, None, 0, "", error=str(e)))
                if show_progress:
                    print(f"[structalign] ({i}/{len(mobiles)}) {s.label}: SKIPPED ({e})", flush=True)
                continue
            out_pdb = out_dir / f"{s.label}_aligned.pdb"
            cmd.save(str(out_pdb), s.label)
            align_result = AlignmentResult(
                s.label, reference_label, result["RMSD"], int(result["alignment_length"]), str(out_pdb),
            )
            results.append(align_result)
            if show_progress:
                print(
                    f"[structalign] ({i}/{len(mobiles)}) {s.label}: rmsd={align_result.rmsd:.3f} "
                    f"({align_result.n_atoms} atoms) -> {out_pdb.name}", flush=True,
                )

    return results
