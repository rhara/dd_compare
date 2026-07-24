"""Fast, offline unit tests for pocketmap.py's cross-protein pocket-residue
classification (no fpocket/PyMOL -- `PocketSelection` is built by hand)."""
from dd_idea.pocket import PocketSelection, Residue
from dd_idea.pocketmap import compare_pocket
from dd_idea.sequence import align_to_reference

REF_SEQ = "ACDEFGHIKLMNPQRSTVWY"


def _pocket(resnums):
    residues = [Residue(chain="A", resnum=n) for n in resnums]
    return PocketSelection(
        receptor_pdb="ref.pdb", fpocket_id=1, rank=1, score=0.5, druggability_score=0.8,
        n_alpha_spheres=10, volume=100.0, residues=residues, center=(0.0, 0.0, 0.0),
        box_center=[0.0, 0.0, 0.0], box_size=[10.0, 10.0, 10.0],
    )


def test_compare_pocket_classifies_identical_conservative_nonconservative_and_gap():
    # Build a target that: keeps position 1 (A) identical, substitutes
    # position 5 (F->D, non-conservative), substitutes position 15 (R->K,
    # conservative), and deletes position 11 (M) entirely (a real
    # insertion/deletion between the two proteins, not a substitution).
    target_list = list(REF_SEQ)
    target_list[4] = "D"  # position 5: F -> D
    target_list[14] = "K"  # position 15: R -> K
    del target_list[10]  # position 11 (M) deleted
    target_seq = "".join(target_list)

    mapping = align_to_reference(REF_SEQ, "REF", target_seq, "TGT")
    pocket = _pocket([1, 5, 11, 15])

    comparisons = {c.reference_position: c for c in compare_pocket(pocket, mapping)}

    assert comparisons[1].conservation == "identical"
    assert comparisons[1].target_residue == "A"

    assert comparisons[5].conservation == "non-conservative"
    assert comparisons[5].target_residue == "D"

    assert comparisons[11].conservation == "gap"
    assert comparisons[11].target_residue is None
    assert comparisons[11].target_position is None

    assert comparisons[15].conservation == "conservative"
    assert comparisons[15].target_residue == "K"
    assert comparisons[15].target_position == 14  # shifted back by 1 after the position-11 deletion


def test_compare_pocket_self_alignment_is_all_identical():
    mapping = align_to_reference(REF_SEQ, "REF", REF_SEQ, "REF")
    pocket = _pocket([1, 7, 20])
    comparisons = compare_pocket(pocket, mapping)
    assert all(c.conservation == "identical" for c in comparisons)
