"""Pfam/InterPro family-membership search over UniProt's REST API --
finding *other* proteins cross-referenced to the same family id as a given
entry (see `similarity.py`, which drives this to propose candidate similar
proteins). Separate from `uniprot/entry.py`'s single-accession fetch.
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from typing import Optional

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"


def family_cross_references(entry: dict) -> list:
    """Every (database, id, entry_name) triple from the Pfam/InterPro
    cross-references UniProt attaches to this entry -- candidate "family"
    identifiers this module searches for other members of."""
    out = []
    for xref in entry.get("uniProtKBCrossReferences", []):
        if xref.get("database") not in ("Pfam", "InterPro"):
            continue
        name = next((p["value"] for p in xref.get("properties", []) if p["key"] == "EntryName"), "")
        out.append((xref["database"], xref["id"], name))
    return out


def count_family_members(database: str, family_id: str, *, taxon_id: Optional[int] = None) -> int:
    """Number of reviewed (Swiss-Prot) UniProt entries cross-referenced to
    this Pfam/InterPro family id, optionally restricted to one organism.
    Uses the `X-Total-Results` response header (a `size=0` query still
    reports the true total there) rather than fetching any rows."""
    query = f"(xref:{database.lower()}-{family_id}) AND (reviewed:true)"
    if taxon_id is not None:
        query += f" AND (organism_id:{taxon_id})"
    url = f"{UNIPROT_SEARCH}?query={urllib.parse.quote(query)}&size=0&fields=accession"
    req = urllib.request.Request(url, headers={"Accept": "application/json"})
    with urllib.request.urlopen(req) as fh:
        total = fh.headers.get("X-Total-Results")
    return int(total) if total is not None else 0


def list_family_members(
    database: str, family_id: str, *, taxon_id: Optional[int] = None, limit: int = 500,
) -> list:
    """Accessions of every reviewed UniProt entry cross-referenced to this
    Pfam/InterPro family id (optionally restricted to one organism), up to
    `limit`. UniProt's search API paginates at 500 rows/request; `limit`
    caps how many pages are followed rather than always draining the
    entire result set, since some family ids resolve to thousands of
    members (see `similarity.py`'s family-selection heuristic, which
    prefers small families specifically to avoid needing this)."""
    query = f"(xref:{database.lower()}-{family_id}) AND (reviewed:true)"
    if taxon_id is not None:
        query += f" AND (organism_id:{taxon_id})"
    accessions: list = []
    url = f"{UNIPROT_SEARCH}?query={urllib.parse.quote(query)}&size={min(500, limit)}&fields=accession"
    while url and len(accessions) < limit:
        req = urllib.request.Request(url, headers={"Accept": "application/json"})
        with urllib.request.urlopen(req) as fh:
            body = json.load(fh)
            link = fh.headers.get("Link")
        accessions.extend(r["primaryAccession"] for r in body.get("results", []))
        url = _next_link(link)
    return accessions[:limit]


def _next_link(link_header: Optional[str]) -> Optional[str]:
    """Parse the RFC 5988 `Link` header UniProt's search API uses for
    pagination and return the `rel="next"` URL, or None on the last page."""
    if not link_header:
        return None
    for part in link_header.split(","):
        if 'rel="next"' in part:
            return part.split(";")[0].strip().strip("<>")
    return None
