"""AlphaFold DB model fetching. Vendored from `dd_seqalign.fetch` (itself
vendored from `dd_prep.fetch`) -- dd_idea only needs these two small
stdlib-only functions, not the RCSB/PDB-entry lookups `dd_seqalign.fetch`
also has (dd_idea's representative structure is always the AlphaFold
model -- see README "Why always AlphaFold, never a PDB structure").
"""
from __future__ import annotations

import json
import urllib.request
from pathlib import Path

AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"


def resolve_afdb_pdb_url(uniprot_id: str) -> str:
    """Look up the current model version's PDB URL for a UniProt accession
    via the AlphaFold DB REST API (model version numbers change over time,
    e.g. v4 -> v6, so this must not be hardcoded)."""
    with urllib.request.urlopen(AFDB_API.format(uniprot_id=uniprot_id.upper())) as fh:
        entries = json.load(fh)
    if not entries:
        raise ValueError(f"AlphaFold DB has no entry for {uniprot_id!r}")
    return entries[0]["pdbUrl"]


def download_afdb(uniprot_id: str, dest: Path) -> str:
    """Fetch the current AlphaFold DB model for a UniProt accession and
    return its text contents. Cached: skipped if `dest` already exists."""
    dest = Path(dest)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        url = resolve_afdb_pdb_url(uniprot_id)
        urllib.request.urlretrieve(url, dest)
    return dest.read_text()
