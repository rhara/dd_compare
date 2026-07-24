"""Fast, offline unit tests for rcsb/fetch.py's
`list_all_structures_at_resolution`, with RCSB network calls mocked out.
Patches target `dd_idea.rcsb.fetch.*` directly (unlike `rcsb/select.py`,
this module calls its own `list_pdb_ids_for_uniprot`/
`fetch_entry_metadata`/`download_pdb`, not imported copies)."""
from unittest.mock import patch

from dd_idea.rcsb import EntryMetadata
from dd_idea.rcsb.fetch import count_structures_for_uniprot, list_all_structures_at_resolution


def _fake_download(pdb_id, dest):
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text("ATOM\n")
    return "ATOM\n"


def test_list_all_structures_at_resolution_downloads_flat_by_default(tmp_path):
    metadata = {"1AAA": EntryMetadata(pdb_id="1AAA", method="X-ray", resolution=1.5, title="")}
    with patch("dd_idea.rcsb.fetch.list_pdb_ids_for_uniprot", return_value=["1AAA"]), \
         patch("dd_idea.rcsb.fetch.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.fetch.download_pdb", side_effect=_fake_download):
        kept = list_all_structures_at_resolution("P00000", tmp_path, show_progress=False)

    assert len(kept) == 1
    assert kept[0]["pdb_path"] == str(tmp_path / "raw_pdb" / "1AAA.pdb")


def test_list_all_structures_at_resolution_downloads_into_subdir_when_given(tmp_path):
    metadata = {"1AAA": EntryMetadata(pdb_id="1AAA", method="X-ray", resolution=1.5, title="")}
    with patch("dd_idea.rcsb.fetch.list_pdb_ids_for_uniprot", return_value=["1AAA"]), \
         patch("dd_idea.rcsb.fetch.fetch_entry_metadata", side_effect=lambda pid: metadata[pid]), \
         patch("dd_idea.rcsb.fetch.download_pdb", side_effect=_fake_download):
        kept = list_all_structures_at_resolution("P00000", tmp_path, subdir="CDK2", show_progress=False)

    assert len(kept) == 1
    assert kept[0]["pdb_path"] == str(tmp_path / "raw_pdb" / "CDK2" / "1AAA.pdb")


def test_count_structures_for_uniprot_counts_without_downloading():
    with patch("dd_idea.rcsb.fetch.list_pdb_ids_for_uniprot", return_value=["1AAA", "1BBB", "1CCC"]) as mocked:
        assert count_structures_for_uniprot("P00000") == 3
    mocked.assert_called_once_with("P00000")


def test_count_structures_for_uniprot_zero_when_none():
    with patch("dd_idea.rcsb.fetch.list_pdb_ids_for_uniprot", return_value=[]):
        assert count_structures_for_uniprot("P00000") == 0
