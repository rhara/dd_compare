"""Fast, offline unit tests for sequence.py's cross-protein alignment logic
(no network, no PyMOL/fpocket)."""
from dd_compare.sequence import align_to_reference, conservation, percent_identity

REF = "ACDEFGHIKLMNPQRSTVWY"  # the 20 standard amino acids, each once, in a fixed order


def test_identical_sequences_are_100_percent_identical_and_map_1_to_1():
    mapping = align_to_reference(REF, "REF", REF, "REF2")
    assert mapping.pct_identity == 100.0
    for pos in range(1, len(REF) + 1):
        assert mapping.ref_to_target_pos[pos] == pos
        assert mapping.target_residue(pos) == mapping.reference_residue(pos)


def test_target_insertion_shifts_downstream_mapping():
    # Insert 3 extra residues into the target right after reference
    # position 3 -- everything before the insertion should map 1:1,
    # everything after should be shifted by +3 in the target's numbering.
    # "WWW" (not "III"/etc.) deliberately avoids repeating a residue that's
    # already adjacent to the insertion point: inserting a run of the same
    # residue already there makes the alignment ambiguous (any of the
    # repeated residues could equally be "the" original one), which isn't
    # what this test is trying to exercise.
    target = REF[:3] + "WWW" + REF[3:]
    mapping = align_to_reference(REF, "REF", target, "TGT")
    assert mapping.pct_identity == 100.0  # the insertion isn't a mismatch, just unmatched target residues
    assert mapping.ref_to_target_pos[1] == 1
    assert mapping.ref_to_target_pos[3] == 3
    assert mapping.ref_to_target_pos[4] == 7
    assert mapping.ref_to_target_pos[len(REF)] == len(REF) + 3


def test_target_deletion_leaves_a_gap_and_shifts_downstream_mapping():
    # Delete reference position 11 ('M') entirely from the target --
    # that position should map to None (no counterpart), and everything
    # after it should shift back by 1.
    target = REF[:10] + REF[11:]
    mapping = align_to_reference(REF, "REF", target, "TGT")
    assert mapping.ref_to_target_pos[10] == 10
    assert mapping.ref_to_target_pos[11] is None
    assert mapping.target_residue(11) is None
    assert mapping.ref_to_target_pos[12] == 11
    assert mapping.ref_to_target_pos[len(REF)] == len(REF) - 1


def test_conservation_identical():
    assert conservation("A", "A") == "identical"


def test_conservation_conservative_substitution():
    # BLOSUM62 R-K = +2 (both positively charged) and D-E = +2 (both negatively charged)
    assert conservation("R", "K") == "conservative"
    assert conservation("D", "E") == "conservative"


def test_conservation_non_conservative_substitution():
    # BLOSUM62 F-D = -3 (aromatic vs. negatively charged)
    assert conservation("F", "D") == "non-conservative"


def test_conservation_gap_when_either_side_missing():
    assert conservation(None, "A") == "gap"
    assert conservation("A", None) == "gap"


def test_percent_identity_identical_sequences():
    assert percent_identity(REF, REF) == 100.0


def test_percent_identity_unrelated_sequences_is_low():
    assert percent_identity(REF, "W" * len(REF)) < 20.0
