"""UniProt entry/sequence and AlphaFold DB model fetching. The AFDB
download/URL-resolution functions are vendored from `dd_seqalign.fetch`
(itself vendored from `dd_prep.fetch`) -- dd_compare only needs these two
small stdlib-only functions plus the UniProt REST calls this module adds
(entry JSON, for the protein name and family cross-references
`similarity.py` needs), not the RCSB/PDB-entry lookups `dd_seqalign.fetch`
also has (dd_compare's representative structure is always the AlphaFold
model -- see README "Why always AlphaFold, never a PDB structure").
"""
from __future__ import annotations

import json
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_ENTRY = "https://rest.uniprot.org/uniprotkb/{accession}.json?fields=protein_name,gene_names,organism_name,xref_pfam,xref_interpro"
UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
AFDB_API = "https://alphafold.ebi.ac.uk/api/prediction/{uniprot_id}"


def fetch_uniprot_fasta(accession: str) -> str:
    """The canonical (isoform 1) amino-acid sequence for a UniProt
    accession, as one contiguous string (FASTA header line stripped)."""
    with urllib.request.urlopen(UNIPROT_FASTA.format(accession=accession.upper())) as fh:
        text = fh.read().decode()
    lines = text.splitlines()
    if not lines or not lines[0].startswith(">"):
        raise ValueError(f"UniProt has no entry for {accession!r}")
    return "".join(lines[1:])


def fetch_uniprot_entry(accession: str) -> dict:
    """Protein name, organism, and Pfam/InterPro cross-references for a
    UniProt accession, as the raw REST JSON (see `similarity.py` for how
    the cross-references are used to find candidate similar proteins)."""
    with urllib.request.urlopen(UNIPROT_ENTRY.format(accession=accession.upper())) as fh:
        return json.load(fh)


def protein_name(entry: dict) -> str:
    """Best-effort human-readable name from a `fetch_uniprot_entry` result."""
    try:
        return entry["proteinDescription"]["recommendedName"]["fullName"]["value"]
    except KeyError:
        pass
    try:
        return entry["proteinDescription"]["submissionNames"][0]["fullName"]["value"]
    except (KeyError, IndexError):
        return entry.get("primaryAccession", "?")


def organism_taxon_id(entry: dict) -> Optional[int]:
    return entry.get("organism", {}).get("taxonId")


def family_cross_references(entry: dict) -> list:
    """Every (database, id, entry_name) triple from the Pfam/InterPro
    cross-references UniProt attaches to this entry -- candidate "family"
    identifiers `similarity.py` searches for other members of."""
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
