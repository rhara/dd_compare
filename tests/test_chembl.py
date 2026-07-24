"""Fast, offline unit tests for chembl/resolve.py's target resolution and
chembl/activity.py's activity counting -- monkeypatches
`urllib.request.urlopen` (no real ChEMBL API access). Patches target
`dd_idea.chembl.resolve.urllib...`/`dd_idea.chembl.activity.urllib...`
(the submodules that actually own the `urllib` import), not
`dd_idea.chembl.urllib...` -- the package's `__init__.py` re-exports the
functions but not the `urllib` reference itself."""
import io
import json

import pytest

from dd_idea import chembl
from dd_idea.chembl import activity as chembl_activity
from dd_idea.chembl import resolve as chembl_resolve

SINGLE_PROTEIN_TARGET = {
    "target_chembl_id": "CHEMBL301",
    "pref_name": "Cyclin-dependent kinase 2",
    "target_type": "SINGLE PROTEIN",
    "target_components": [
        {"accession": "P24941", "component_type": "PROTEIN", "relationship": "SINGLE PROTEIN"},
    ],
}

TARGET_SEARCH_RESULTS = {
    "targets": [
        {"target_chembl_id": "CHEMBL301", "pref_name": "Cyclin-dependent kinase 2", "organism": "Homo sapiens"},
    ],
}

EMPTY_TARGET_SEARCH_RESULTS = {"targets": []}

ACTIVITY_PAGE = {"page_meta": {"limit": 1, "offset": 0, "total_count": 3015, "next": None, "previous": None}, "activities": []}

MULTI_COMPONENT_TARGET = {
    "target_chembl_id": "CHEMBL2094127",
    "pref_name": "Some protein complex",
    "target_type": "PROTEIN COMPLEX",
    "target_components": [
        {"accession": "P11111", "component_type": "PROTEIN", "relationship": "PROTEIN COMPLEX"},
        {"accession": "P22222", "component_type": "PROTEIN", "relationship": "PROTEIN COMPLEX"},
    ],
}


class _FakeResponse:
    def __init__(self, payload: dict):
        self._body = json.dumps(payload).encode()

    def __enter__(self):
        return io.BytesIO(self._body)

    def __exit__(self, *exc):
        return False


def test_resolve_chembl_target_to_uniprot_single_protein(monkeypatch):
    monkeypatch.setattr(chembl_resolve.urllib.request, "urlopen", lambda url: _FakeResponse(SINGLE_PROTEIN_TARGET))
    assert chembl.resolve_chembl_target_to_uniprot("CHEMBL301") == "P24941"


def test_resolve_chembl_target_to_uniprot_lowercases_are_accepted(monkeypatch):
    monkeypatch.setattr(chembl_resolve.urllib.request, "urlopen", lambda url: _FakeResponse(SINGLE_PROTEIN_TARGET))
    assert chembl.resolve_chembl_target_to_uniprot("chembl301") == "P24941"


def test_resolve_chembl_target_to_uniprot_rejects_multi_component_targets(monkeypatch):
    monkeypatch.setattr(chembl_resolve.urllib.request, "urlopen", lambda url: _FakeResponse(MULTI_COMPONENT_TARGET))
    with pytest.raises(ValueError):
        chembl.resolve_chembl_target_to_uniprot("CHEMBL2094127")


def test_resolve_uniprot_to_chembl_targets_returns_matches(monkeypatch):
    monkeypatch.setattr(chembl_resolve.urllib.request, "urlopen", lambda url: _FakeResponse(TARGET_SEARCH_RESULTS))
    targets = chembl.resolve_uniprot_to_chembl_targets("P24941")
    assert targets == [{"target_chembl_id": "CHEMBL301", "pref_name": "Cyclin-dependent kinase 2", "organism": "Homo sapiens"}]


def test_resolve_uniprot_to_chembl_targets_empty_when_no_target(monkeypatch):
    monkeypatch.setattr(chembl_resolve.urllib.request, "urlopen", lambda url: _FakeResponse(EMPTY_TARGET_SEARCH_RESULTS))
    assert chembl.resolve_uniprot_to_chembl_targets("Q8IZL9") == []


def test_count_activities_reads_total_count_without_pagination(monkeypatch):
    monkeypatch.setattr(chembl_activity.urllib.request, "urlopen", lambda url: _FakeResponse(ACTIVITY_PAGE))
    assert chembl.count_activities("CHEMBL301") == 3015


def test_count_activities_default_assay_type_is_binding_only():
    # DEFAULT_ASSAY_TYPES is what count_activities uses unless overridden --
    # matches dd_chembl.fetch.DEFAULT_ASSAY_TYPES exactly, so the count
    # reflects what dd_chembl would actually pull for QSAR training.
    assert chembl.DEFAULT_ASSAY_TYPES == ("B",)
