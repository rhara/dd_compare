"""Single-UniProt-entry fetching (sequence + entry JSON) and accessors for
the fields dd_idea's tables/reports actually read out of that JSON. See
`uniprot/family.py` for the separate Pfam/InterPro-family-membership
search UniProt's REST API also offers.
"""
from __future__ import annotations

import json
import urllib.request
from typing import Optional

UNIPROT_FASTA = "https://rest.uniprot.org/uniprotkb/{accession}.fasta"
UNIPROT_ENTRY = (
    "https://rest.uniprot.org/uniprotkb/{accession}.json"
    "?fields=protein_name,gene_names,organism_name,length,cc_similarity,xref_pfam,xref_interpro"
)


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
    """Protein name, organism, length, family (SIMILARITY comment), and
    Pfam/InterPro cross-references for a UniProt accession, as the raw
    REST JSON (see `uniprot/family.py` for how the cross-references are
    used to find candidate similar proteins)."""
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


def organism_name(entry: dict) -> str:
    return entry.get("organism", {}).get("scientificName", "-")


def gene_name(entry: dict) -> str:
    genes = entry.get("genes", [])
    if genes and "geneName" in genes[0]:
        return genes[0]["geneName"]["value"]
    return "-"


def canonical_length(entry: dict) -> Optional[int]:
    return entry.get("sequence", {}).get("length")


def family_string(entry: dict) -> str:
    """The UniProt SIMILARITY comment -- family/subfamily classification
    in UniProt's own words, e.g. "Belongs to the protein kinase
    superfamily. CMGC Ser/Thr protein kinase family. CDC2/CDKX subfamily"
    -- or `"(none reported)"` if this entry has no such comment. This is
    prose, not a structured id like `family_cross_references`'
    Pfam/InterPro ids, but it's the single most legible "is this actually
    related" signal when scanning a table of BLAST hits."""
    for comment in entry.get("comments", []):
        if comment.get("commentType") == "SIMILARITY":
            texts = comment.get("texts", [])
            if texts:
                return texts[0]["value"]
    return "(none reported)"
