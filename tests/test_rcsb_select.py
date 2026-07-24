"""Fast, offline unit tests for rcsb/select.py's `select_pdb_structures`
ligand-diversity/resolution-fallback decision logic and
`rehydrate_selection`, with RCSB network calls mocked out. Patches target
`dd_idea.rcsb.select.*` (the names as `select.py` imported them from
`rcsb/fetch.py`), not `dd_idea.rcsb.fetch.*`, since that's where
`select_pdb_structures` actually looks them up."""
from unittest.mock import patch

from dd_idea.rcsb import EntryMetadata, rehydrate_selection, select_pdb_structures

_APO_PDB = "ATOM      1  CA  ALA A   1      11.104   6.134  -6.504  1.00 20.00           C\n"


def _hetatm_line(serial: int, atom_name: str, resname: str, resseq: int, x: float, y: float, z: float) -> str:
    # Fixed-column PDB format -- placed at exact 0-indexed offsets (matching
    # pdbio.py's own slicing: resname at [17:20], chain at [21], resseq at
    # [22:26], coords at [30:38]/[38:46]/[46:54]) rather than composed via
    # f-string field widths, which are easy to miscount by one column and
    # silently shift every field after the mistake.
    line = list(" " * 80)
    line[0:6] = "HETATM"
    line[6:11] = f"{serial:>5}"
    line[12:16] = f"{atom_name:<4}"
    line[17:20] = f"{resname:>3}"
    line[21] = "A"
    line[22:26] = f"{resseq:>4}"
    line[30:38] = f"{x:>8.3f}"
    line[38:46] = f"{y:>8.3f}"
    line[46:54] = f"{z:>8.3f}"
    line[76:78] = " C"
    return "".join(line) + "\n"


def _ligand_pdb(resname: str) -> str:
    assert len(resname) == 3, "PDB resnames are fixed-width 3 characters"
    return _APO_PDB + "".join(
        _hetatm_line(i, f"C{i}", resname, 201, 10.0 + i, 7.0 + i, -5.0 + i) for i in range(2, 7)
    )


def test_select_pdb_structures_prefers_ligand_bound_over_better_resolution(tmp_path):
    # 2ABC (apo, best resolution) should lose to 3XYZ (worse resolution but
    # ligand-bound) -- resolution alone is not the deciding factor.
    ids = ["2ABC", "3XYZ"]
    metadata = {
        "2ABC": EntryMetadata(pdb_id="2ABC", method="X-ray", resolution=1.0, title=""),
        "3XYZ": EntryMetadata(pdb_id="3XYZ", method="X-ray", resolution=1.8, title=""),
    }
    pdb_text = {"2ABC": _APO_PDB, "3XYZ": _ligand_pdb("LIG")}

    def fake_download(pdb_id, dest):
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(pdb_text[pdb_id])
        return pdb_text[pdb_id]

    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_idea.rcsb.select.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.select.download_pdb", side_effect=fake_download):
        sels = select_pdb_structures("P00000", "A", tmp_path, scan_cap=10, show_progress=False)

    assert [s.pdb_id for s in sels] == ["3XYZ"]
    assert sels[0].ligand_resname == "LIG"


def test_select_pdb_structures_falls_back_to_best_resolution_when_all_apo(tmp_path):
    ids = ["2ABC", "3XYZ"]
    metadata = {
        "2ABC": EntryMetadata(pdb_id="2ABC", method="X-ray", resolution=1.0, title=""),
        "3XYZ": EntryMetadata(pdb_id="3XYZ", method="X-ray", resolution=1.8, title=""),
    }

    def fake_download(pdb_id, dest):
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(_APO_PDB)
        return _APO_PDB

    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_idea.rcsb.select.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.select.download_pdb", side_effect=fake_download):
        sels = select_pdb_structures("P00000", "A", tmp_path, scan_cap=10, show_progress=False)

    assert len(sels) == 1
    assert sels[0].pdb_id == "2ABC"  # best resolution among the (all-apo) candidates
    assert sels[0].ligand_resname is None


def test_select_pdb_structures_excludes_worse_than_resolution_cutoff(tmp_path):
    # 1AAA is ligand-bound but resolution 3.0 -- worse than the default
    # 2.0Å cutoff, so it must be skipped entirely (never even downloaded);
    # 2BBB (apo, 1.5Å) is the only one that meets the bar, so it wins the
    # apo fallback despite being scanned second.
    ids = ["1AAA", "2BBB"]
    metadata = {
        "1AAA": EntryMetadata(pdb_id="1AAA", method="X-ray", resolution=3.0, title=""),
        "2BBB": EntryMetadata(pdb_id="2BBB", method="X-ray", resolution=1.5, title=""),
    }
    pdb_text = {"1AAA": _ligand_pdb("LGA"), "2BBB": _APO_PDB}
    downloaded = []

    def fake_download(pdb_id, dest):
        downloaded.append(pdb_id)
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(pdb_text[pdb_id])
        return pdb_text[pdb_id]

    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_idea.rcsb.select.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.select.download_pdb", side_effect=fake_download):
        sels = select_pdb_structures("P00000", "A", tmp_path, scan_cap=10, show_progress=False)

    assert [s.pdb_id for s in sels] == ["2BBB"]
    assert "1AAA" not in downloaded  # excluded by the resolution cutoff before ever being downloaded


def test_select_pdb_structures_excludes_unreported_resolution(tmp_path):
    # An NMR-style entry with no reported resolution at all must be
    # excluded regardless of how loose the cutoff is (None can't satisfy
    # a "<= cutoff" comparison).
    ids = ["1NMR"]
    metadata = {"1NMR": EntryMetadata(pdb_id="1NMR", method="NMR", resolution=None, title="")}

    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_idea.rcsb.select.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]):
        sels = select_pdb_structures(
            "P00000", "A", tmp_path, scan_cap=10, resolution_cutoff=100.0, show_progress=False,
        )

    assert sels == []


def test_select_pdb_structures_resolution_cutoff_is_configurable(tmp_path):
    # The same 3.0Å ligand-bound entry that the default cutoff excludes
    # (see test above) is accepted once the caller loosens the cutoff.
    ids = ["1AAA"]
    metadata = {"1AAA": EntryMetadata(pdb_id="1AAA", method="X-ray", resolution=3.0, title="")}
    pdb_text = {"1AAA": _ligand_pdb("LGA")}

    def fake_download(pdb_id, dest):
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(pdb_text[pdb_id])
        return pdb_text[pdb_id]

    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_idea.rcsb.select.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.select.download_pdb", side_effect=fake_download):
        sels = select_pdb_structures(
            "P00000", "A", tmp_path, scan_cap=10, resolution_cutoff=3.5, show_progress=False,
        )

    assert [s.pdb_id for s in sels] == ["1AAA"]


def test_select_pdb_structures_returns_empty_when_no_structures_exist(tmp_path):
    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=[]):
        sels = select_pdb_structures("P00000", "A", tmp_path, show_progress=False)
    assert sels == []


def test_select_pdb_structures_collects_multiple_distinct_ligands_up_to_cap(tmp_path):
    # Five ligand-bound candidates, best-resolution-first: LGA, LGB, LGA
    # again (a duplicate ligand from a worse-resolution entry -- should be
    # skipped in favor of chemical diversity, not counted toward the cap),
    # then LGC, then LGD. With max_structures=3, the 3 distinct ligands
    # (LGA, LGB, LGC) are kept, 3CCC's duplicate LGA is skipped, and 5EEE
    # is never even scanned once the cap is reached.
    ids = ["1AAA", "2BBB", "3CCC", "4DDD", "5EEE"]
    metadata = {pid: EntryMetadata(pdb_id=pid, method="X-ray", resolution=float(i), title="")
                for i, pid in enumerate(ids, start=1)}
    pdb_text = {"1AAA": _ligand_pdb("LGA"), "2BBB": _ligand_pdb("LGB"),
                "3CCC": _ligand_pdb("LGA"), "4DDD": _ligand_pdb("LGC"), "5EEE": _ligand_pdb("LGD")}
    downloaded = []

    def fake_download(pdb_id, dest):
        downloaded.append(pdb_id)
        dest = tmp_path / f"{pdb_id}.pdb"
        dest.write_text(pdb_text[pdb_id])
        return pdb_text[pdb_id]

    with patch("dd_idea.rcsb.select.list_pdb_ids_for_uniprot", return_value=ids), \
         patch("dd_idea.rcsb.select.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.select.download_pdb", side_effect=fake_download):
        sels = select_pdb_structures(
            "P00000", "A", tmp_path, scan_cap=10, max_structures=3, resolution_cutoff=10.0, show_progress=False,
        )  # resolution_cutoff loosened -- this test is about dedup/cap logic, not resolution filtering

    assert [s.pdb_id for s in sels] == ["1AAA", "2BBB", "4DDD"]
    assert [s.ligand_resname for s in sels] == ["LGA", "LGB", "LGC"]
    assert "3CCC" in downloaded  # scanned and correctly recognized as a duplicate ligand
    assert "5EEE" not in downloaded  # stopped scanning once max_structures distinct ligands were reached


def test_rehydrate_selection_reconstructs_without_network(tmp_path):
    # rehydrate_selection is what pipeline.analyze uses to reload a
    # fetch-time pick from manifest.json -- no RCSB calls at all, just
    # re-reading the already-downloaded file and re-running the local
    # canonical-sequence alignment.
    pdb_path = tmp_path / "9XYZ.pdb"
    pdb_path.write_text(_ligand_pdb("LIG"))

    sel = rehydrate_selection(
        "P00000", "A", pdb_id="9XYZ", resolution=1.23, chain_id="A",
        pdb_path=str(pdb_path), ligand_resname="LIG",
    )

    assert sel.pdb_id == "9XYZ"
    assert sel.resolution == 1.23
    assert sel.ligand_resname == "LIG"
    assert sel.chain_alignment.resseq_for_canonical(1) == 1  # the fixture's single ALA at resseq 1
