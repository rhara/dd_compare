"""ChEMBL target-ID -> UniProt-accession resolution, for `dd_idea-search`'s
ChEMBL-target-ID input form. Talks to the ChEMBL REST API directly over
`urllib.request` (same convention as `dd_chembl/dd_chembl/fetch.py`,
including its `CHEMBL_API_BASE`), rather than depending on
chembl_webresource_client -- vendored, not imported, per this whole
project family's convention of not depending on any sibling `dd_*` package
at install time (dd_idea only needs this one direction -- ChEMBL target ID
in, UniProt accession out -- not dd_chembl's much larger bioactivity-data
pipeline).
"""
from __future__ import annotations

import json
import urllib.request

CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
CHEMBL_TARGET = CHEMBL_API_BASE + "/target/{chembl_id}.json"


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
