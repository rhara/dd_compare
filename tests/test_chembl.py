"""Fast, offline unit tests for chembl/resolve.py's target resolution --
monkeypatches `urllib.request.urlopen` (no real ChEMBL API access). Patches
target `dd_idea.chembl.resolve.urllib...` (the submodule that actually
owns the `urllib` import), not `dd_idea.chembl.urllib...` -- the package's
`__init__.py` re-exports the function but not the `urllib` reference
itself."""
import io
import json

import pytest

from dd_idea import chembl
from dd_idea.chembl import resolve as chembl_resolve

SINGLE_PROTEIN_TARGET = {
    "target_chembl_id": "CHEMBL301",
    "pref_name": "Cyclin-dependent kinase 2",
    "target_type": "SINGLE PROTEIN",
    "target_components": [
        {"accession": "P24941", "component_type": "PROTEIN", "relationship": "SINGLE PROTEIN"},
    ],
}

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
