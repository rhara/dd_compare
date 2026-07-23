"""Fast, offline unit tests for similarity.py's family-selection and
candidate-ranking logic -- `fetch`'s network calls are all monkeypatched
(no real UniProt access)."""
from dd_compare import fetch, similarity

SEQ = "ACDEFGHIKLMNPQRSTVWY"


def test_choose_family_picks_smallest_family_with_more_than_one_member(monkeypatch):
    monkeypatch.setattr(
        fetch, "family_cross_references",
        lambda entry: [("Pfam", "PF00001", "Broad"), ("InterPro", "IPR00002", "Narrow"), ("InterPro", "IPR00003", "TooNarrow")],
    )
    counts = {("Pfam", "PF00001"): 500, ("InterPro", "IPR00002"): 10, ("InterPro", "IPR00003"): 1}
    monkeypatch.setattr(fetch, "count_family_members", lambda db, fid, taxon_id=None: counts[(db, fid)])

    choice = similarity.choose_family({}, taxon_id=9606, show_progress=False)
    assert (choice.database, choice.family_id, choice.member_count) == ("InterPro", "IPR00002", 10)


def test_choose_family_prefers_interpro_over_pfam_on_a_count_tie(monkeypatch):
    monkeypatch.setattr(
        fetch, "family_cross_references", lambda entry: [("Pfam", "PF00001", "P"), ("InterPro", "IPR00001", "I")],
    )
    monkeypatch.setattr(fetch, "count_family_members", lambda db, fid, taxon_id=None: 5)

    choice = similarity.choose_family({}, taxon_id=9606, show_progress=False)
    assert choice.database == "InterPro"


def test_choose_family_returns_none_when_every_family_is_a_singleton(monkeypatch):
    monkeypatch.setattr(fetch, "family_cross_references", lambda entry: [("Pfam", "PF00001", "Self only")])
    monkeypatch.setattr(fetch, "count_family_members", lambda db, fid, taxon_id=None: 1)

    assert similarity.choose_family({}, taxon_id=9606, show_progress=False) is None


def test_discover_ranks_by_identity_and_excludes_the_seed_itself(monkeypatch):
    monkeypatch.setattr(fetch, "fetch_uniprot_entry", lambda acc: {"accession": acc})
    monkeypatch.setattr(fetch, "protein_name", lambda entry: f"Protein {entry['accession']}")
    monkeypatch.setattr(fetch, "organism_taxon_id", lambda entry: 9606)
    monkeypatch.setattr(fetch, "family_cross_references", lambda entry: [("InterPro", "IPR00001", "Fam")])
    monkeypatch.setattr(fetch, "count_family_members", lambda db, fid, taxon_id=None: 3)
    monkeypatch.setattr(
        fetch, "list_family_members", lambda db, fid, taxon_id=None, limit=500: ["SEED", "CLOSE", "FAR"],
    )
    seqs = {"SEED": SEQ, "CLOSE": SEQ, "FAR": "W" * len(SEQ)}
    monkeypatch.setattr(fetch, "fetch_uniprot_fasta", lambda acc: seqs[acc])

    result = similarity.discover("SEED", max_candidates=5, show_progress=False)

    accessions = [c["accession"] for c in result["candidates"]]
    assert "SEED" not in accessions
    assert accessions[0] == "CLOSE"
    assert result["candidates"][0]["pct_identity"] == 100.0


def test_discover_respects_max_candidates(monkeypatch):
    monkeypatch.setattr(fetch, "fetch_uniprot_entry", lambda acc: {"accession": acc})
    monkeypatch.setattr(fetch, "protein_name", lambda entry: entry["accession"])
    monkeypatch.setattr(fetch, "organism_taxon_id", lambda entry: 9606)
    monkeypatch.setattr(fetch, "family_cross_references", lambda entry: [("InterPro", "IPR00001", "Fam")])
    monkeypatch.setattr(fetch, "count_family_members", lambda db, fid, taxon_id=None: 4)
    monkeypatch.setattr(
        fetch, "list_family_members", lambda db, fid, taxon_id=None, limit=500: ["SEED", "A", "B", "C"],
    )
    monkeypatch.setattr(fetch, "fetch_uniprot_fasta", lambda acc: SEQ)

    result = similarity.discover("SEED", max_candidates=2, show_progress=False)
    assert len(result["candidates"]) == 2
