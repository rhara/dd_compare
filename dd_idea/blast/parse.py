"""Turn a BLASTP-vs-Swiss-Prot XML result (from `blast/query.py`, or any
previously-cached `raw_blast/blastp_swissprot.xml`) into ranked
`BlastHit`s -- pure parsing, no network access, so it can run against an
already-fetched XML with no NCBI round-trip at all."""
from __future__ import annotations

import io
import re
from dataclasses import dataclass
from typing import List, Optional

from Bio.Blast import NCBIXML

# Swiss-Prot BLAST hit ids look like "sp|P24941.1|CDK2_HUMAN"
_SWISSPROT_HIT_ID_RE = re.compile(r"^sp\|([A-Z0-9]+)(?:\.\d+)?\|")
_UNIPROT_ACC_RE = re.compile(r"^[A-Z][A-Z0-9]{5}([A-Z0-9]{4})?$")


@dataclass
class BlastHit:
    accession: str
    description: str
    pct_identity: float
    evalue: float
    align_length: int

    def to_dict(self) -> dict:
        return {
            "accession": self.accession, "description": self.description,
            "pct_identity": self.pct_identity, "evalue": self.evalue, "align_length": self.align_length,
        }


def _extract_accession(alignment) -> Optional[str]:
    """Prefer Biopython's own `alignment.accession` (parsed straight from
    the XML's `<Hit_accession>`, version-suffix-free); fall back to
    regexing `hit_id` (e.g. `sp|P24941.1|CDK2_HUMAN`) if that's missing or
    doesn't look like a real UniProt accession."""
    acc = (getattr(alignment, "accession", "") or "").upper()
    if _UNIPROT_ACC_RE.match(acc):
        return acc
    m = _SWISSPROT_HIT_ID_RE.match(alignment.hit_id)
    return m.group(1).upper() if m else None


def parse_blast_hits(xml_text: str, *, exclude_accession: Optional[str] = None, max_hits: int = 100) -> List[BlastHit]:
    """One `BlastHit` per Swiss-Prot alignment (best HSP only), ranked by
    %identity -- NCBI's own hit order is already e-value-ranked, but ties
    at very low e-values are common among close paralogs."""
    blast_record = NCBIXML.read(io.StringIO(xml_text))
    hits: List[BlastHit] = []
    seen: set = set()
    exclude = exclude_accession.upper() if exclude_accession else None
    for alignment in blast_record.alignments:
        acc = _extract_accession(alignment)
        if acc is None:
            continue  # not a Swiss-Prot hit with a recognizable accession -- skip defensively
        if acc == exclude or acc in seen:
            continue
        hsp = alignment.hsps[0]  # best-scoring HSP for this alignment
        pct_identity = 100.0 * hsp.identities / hsp.align_length
        seen.add(acc)
        hits.append(BlastHit(
            accession=acc, description=alignment.hit_def, pct_identity=pct_identity,
            evalue=hsp.expect, align_length=hsp.align_length,
        ))
    hits.sort(key=lambda h: h.pct_identity, reverse=True)
    return hits[:max_hits]
