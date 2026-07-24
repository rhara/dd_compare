"""ChEMBL target <-> UniProt-accession resolution, both directions. Talks
to the ChEMBL REST API directly over `urllib.request` (same convention as
`dd_chembl/dd_chembl/fetch.py`, including its `CHEMBL_API_BASE`), rather
than depending on chembl_webresource_client -- vendored, not imported, per
this whole project family's convention of not depending on any sibling
`dd_*` package at install time (dd_idea only needs target resolution, not
dd_chembl's much larger bioactivity-data pipeline; see `chembl/activity.py`
for the one piece of that dd_idea does need -- activity counts).
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import List

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_TARGET = CHEMBL_API_BASE + "/target/{chembl_id}.json"
CHEMBL_TARGET_SEARCH = CHEMBL_API_BASE + "/target.json"


def resolve_chembl_target_to_uniprot(chembl_id: str) -> str:
    """The single UniProt accession behind a ChEMBL target ID. Only
    `target_type == "SINGLE PROTEIN"` targets with exactly one `PROTEIN`
    target_component are supported -- multi-component targets (protein
    complexes, protein families) have no single canonical sequence for
    the rest of this pipeline (BLAST/AlphaFold/pocket detection) to anchor
    on, so they're rejected with a clear error rather than guessing which
    component the caller meant."""
    chembl_id = chembl_id.upper()
    with urllib.request.urlopen(CHEMBL_TARGET.format(chembl_id=chembl_id)) as fh:
        target = json.load(fh)

    target_type = target.get("target_type")
    components = [c for c in target.get("target_components", []) if c.get("component_type") == "PROTEIN"]
    if target_type != "SINGLE PROTEIN" or len(components) != 1:
        pref_name = target.get("pref_name", "?")
        raise ValueError(
            f"{chembl_id} ({pref_name!r}, target_type={target_type!r}, {len(components)} protein component(s)) "
            f"isn't a single-protein target -- dd_idea-search needs one canonical UniProt accession to BLAST "
            f"from. Pass the specific UniProt accession you want directly instead."
        )
    accession = components[0].get("accession")
    if not accession:
        raise ValueError(f"{chembl_id}: its one protein component has no UniProt accession on record")
    return accession


def resolve_uniprot_to_chembl_targets(accession: str) -> List[dict]:
    """Every ChEMBL `SINGLE PROTEIN` target for a UniProt accession (a
    UniProt accession occasionally maps to more than one -- e.g. distinct
    isoform/mutant target definitions -- so all are returned), as
    `{target_chembl_id, pref_name, organism}` dicts. Same query dd_chembl's
    own `resolve_targets` uses; `[]` if ChEMBL has no single-protein target
    for this accession at all."""
    query = urllib.parse.urlencode({
        "target_components__accession": accession.upper(),
        "target_type": "SINGLE PROTEIN",
        "only": "target_chembl_id,pref_name,organism",
        "format": "json",
    })
    with urllib.request.urlopen(f"{CHEMBL_TARGET_SEARCH}?{query}") as fh:
        body = json.load(fh)
    return [
        {"target_chembl_id": t["target_chembl_id"], "pref_name": t.get("pref_name", "?"), "organism": t.get("organism", "?")}
        for t in body.get("targets", [])
    ]
