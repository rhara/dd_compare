"""Fast, offline unit tests for pdbstruct.py.

The chain-alignment/picking tests are ported alongside the vendored code
they cover (from dd_seqalign.sequence's own test_sequence.py -- see
pdbstruct.py's module docstring), following this project family's standing
pattern of porting a vendored module's test coverage, not just its code.
`select_pdb_structure`'s ligand-vs-resolution-fallback decision is new to
dd_compare, tested here with the network calls mocked out."""
from unittest.mock import patch

from dd_compare.pdbstruct import (
    ChainSequence,
    EntryMetadata,
    align_to_canonical,
    pick_target_chain,
    select_pdb_structure,
)

CANONICAL = "MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLDSI"


def test_align_to_canonical_full_coverage_no_mismatch():
    chain = ChainSequence(chain_id="A", residues=list(enumerate(CANONICAL, start=1)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_covered == len(CANONICAL)
    assert aln.n_mismatch == 0
    assert aln.resseq_for_canonical(1) == 1


def test_align_to_canonical_fragment_has_free_end_gaps():
    fragment = CANONICAL[10:40]  # a co-crystal that only resolved residues 11-40
    chain = ChainSequence(chain_id="A", residues=list(enumerate(fragment, start=11)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_covered == 30
    assert aln.resseq_for_canonical(1) is None  # outside the fragment: missing, not a bad alignment
    assert aln.resseq_for_canonical(11) == 11
    assert aln.resseq_for_canonical(40) == 40


def test_align_to_canonical_flags_point_mutation_not_missing():
    mutant = CANONICAL[:20] + "A" + CANONICAL[21:]
    chain = ChainSequence(chain_id="A", residues=list(enumerate(mutant, start=1)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_mismatch == 1
    mismatch = next(r for r in aln.residues if r.status == "mismatch")
    assert mismatch.canonical_pos == 21


def test_pick_target_chain_prefers_identity_over_raw_coverage():
    # A bound partner chain (e.g. a cyclin, a CAK) can align across the same
    # span as the true target chain while being mostly mismatches --
    # pick_target_chain must rank by matching residues, not raw coverage.
    true_chain = ChainSequence(chain_id="A", residues=list(enumerate(CANONICAL, start=1)))
    partner_seq = "".join("X" if i % 2 == 0 else c for i, c in enumerate(CANONICAL))
    partner_chain = ChainSequence(chain_id="B", residues=list(enumerate(partner_seq, start=1)))

    alignments = {
        "A": align_to_canonical(true_chain, CANONICAL),
        "B": align_to_canonical(partner_chain, CANONICAL),
    }
    assert pick_target_chain(alignments) == "A"


_APO_PDB = "ATOM      1  CA  ALA A   1      11.104   6.134  -6.504  1.00 20.00           C\n"
_LIGAND_PDB = _APO_PDB + "HETATM    2  C1  LIG A 201      12.000   7.000  -5.000  1.00 20.00           C\n" \
    "HETATM    3  C2  LIG A 201      13.000   8.000  -4.000  1.00 20.00           C\n" \
    "HETATM    4  C3  LIG A 201      14.000   9.000  -3.000  1.00 20.00           C\n" \
    "HETATM    5  C4  LIG A 201      15.000  10.000  -2.000  1.00 20.00           C\n" \
    "HETATM    6  C5  LIG A 201      16.000  11.000  -1.000  1.00 20.00           C\n"


def test_select_pdb_structure_prefers_ligand_bound_over_better_resolution(tmp_path):
    # 2ABC (apo, best resolution) should lose to 3XYZ (worse resolution but
    # ligand-bound) -- resolution alone is not the deciding factor.
    ids = ["2ABC", "3XYZ"]
    metadata = {
        "2ABC": EntryMetadata(pdb_id="2ABC", method="X-ray", resolution=1.0, title=""),
        "3XYZ": EntryMetadata(pdb_id="3XYZ", method="X-ray", resolution=2.5, title=""),
    }
    pdb_text = {"2ABC": _APO_PDB, "3XYZ": _LIGAND_PDB}

    def fake_download(pdb_id, dest):
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(pdb_text[pdb_id])
        return pdb_text[pdb_id]

    with patch("dd_compare.pdbstruct.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_compare.pdbstruct.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_compare.pdbstruct.download_pdb", side_effect=fake_download):
        sel = select_pdb_structure("P00000", "A", tmp_path, scan_cap=10, show_progress=False)

    assert sel.pdb_id == "3XYZ"
    assert sel.ligand_resname == "LIG"


def test_select_pdb_structure_falls_back_to_best_resolution_when_all_apo(tmp_path):
    ids = ["2ABC", "3XYZ"]
    metadata = {
        "2ABC": EntryMetadata(pdb_id="2ABC", method="X-ray", resolution=1.0, title=""),
        "3XYZ": EntryMetadata(pdb_id="3XYZ", method="X-ray", resolution=2.5, title=""),
    }

    def fake_download(pdb_id, dest):
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(_APO_PDB)
        return _APO_PDB

    with patch("dd_compare.pdbstruct.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_compare.pdbstruct.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_compare.pdbstruct.download_pdb", side_effect=fake_download):
        sel = select_pdb_structure("P00000", "A", tmp_path, scan_cap=10, show_progress=False)

    assert sel.pdb_id == "2ABC"  # best resolution among the (all-apo) candidates
    assert sel.ligand_resname is None


def test_select_pdb_structure_returns_none_when_no_structures_exist(tmp_path):
    with patch("dd_compare.pdbstruct.list_pdb_ids_for_uniprot", return_value=[]):
        sel = select_pdb_structure("P00000", "A", tmp_path, show_progress=False)
    assert sel is None
