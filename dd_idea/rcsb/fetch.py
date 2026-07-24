"""RCSB PDB lookup/download: which real structures (if any) exist for a
UniProt accession, their metadata, and downloading them. Pure RCSB API
access -- no chain alignment, no selection heuristics (see `rcsb/chain.py`
and `rcsb/select.py` for those).
"""
from __future__ import annotations

import json
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

RCSB_SEARCH = "https://search.rcsb.org/rcsbsearch/v2/query"
RCSB_ENTRY = "https://data.rcsb.org/rest/v1/core/entry/{pdb_id}"
RCSB_PDB = "https://files.rcsb.org/download/{pdb_id}.pdb"


def list_pdb_ids_for_uniprot(accession: str) -> List[str]:
    """Every RCSB PDB entry ID cross-referenced (via SIFTS) to this UniProt
    accession, best-resolution-first (server-side sort) -- a well-studied
    target can have hundreds of entries (e.g. CDK2's 512), and sorting on
    RCSB's side means callers only need to fetch per-entry metadata/
    coordinates for as many of them as they actually scan, instead of
    every one just to rank them."""
    query = {
        "query": {
            "type": "terminal",
            "service": "text",
            "parameters": {
                "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
                "operator": "exact_match",
                "value": accession.upper(),
            },
        },
        "return_type": "entry",
        "request_options": {
            "return_all_hits": True,
            "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}],
        },
    }
    req = urllib.request.Request(
        RCSB_SEARCH, data=json.dumps(query).encode(), headers={"Content-Type": "application/json"}
    )
    with urllib.request.urlopen(req) as fh:
        body = fh.read()
    if not body:
        return []  # RCSB responds 204 No Content (empty body) when a UniProt accession has zero structures
    return [hit["identifier"] for hit in json.loads(body).get("result_set", [])]


@dataclass
class EntryMetadata:
    pdb_id: str
    method: str
    resolution: Optional[float]
    title: str


def fetch_entry_metadata(pdb_id: str) -> EntryMetadata:
    """Experimental method and resolution (None for methods that don't
    report one, e.g. NMR), straight from RCSB's entry-level summary -- cheap
    (no coordinates), used to rank candidates before downloading any of
    them."""
    with urllib.request.urlopen(RCSB_ENTRY.format(pdb_id=pdb_id.upper())) as fh:
        entry = json.load(fh)
    info = entry.get("rcsb_entry_info", {})
    resolution_list = info.get("resolution_combined") or []
    return EntryMetadata(
        pdb_id=pdb_id.upper(),
        method=info.get("experimental_method", "unknown"),
        resolution=resolution_list[0] if resolution_list else None,
        title=entry.get("struct", {}).get("title", ""),
    )


def download_pdb(pdb_id: str, dest: Path) -> str:
    """Fetch a raw PDB entry from RCSB and return its text contents.
    Cached: skipped if `dest` already exists."""
    dest = Path(dest)
    if not dest.exists():
        dest.parent.mkdir(parents=True, exist_ok=True)
        urllib.request.urlretrieve(RCSB_PDB.format(pdb_id=pdb_id.upper()), dest)
    return dest.read_text()


def list_all_structures_at_resolution(
    accession: str, out_dir: Union[str, Path], *, resolution_cutoff: float = 2.0, show_progress: bool = True,
) -> List[dict]:
    """Every RCSB structure of `accession` at or better than
    `resolution_cutoff` -- no cap, no ligand preference, no dedup by
    ligand (unlike `rcsb.select.select_pdb_structures`, which picks up to
    a handful of distinct-ligand entries for a single illustrative
    overlay). Used by `search.py` to gather every reasonable-quality
    template for pocket detection/restrained-MD, where more structural
    diversity is better, not worse. Returns plain dicts (`pdb_id`/
    `resolution`/`method`/`title`/`pdb_path`), not `SelectedPdbStructure`
    -- no canonical-sequence chain alignment is done here, since bulk
    template gathering doesn't need per-residue numbering the way the
    reference-pocket overlay does."""
    out_dir = Path(out_dir)
    raw_pdb_dir = out_dir / "raw_pdb"
    raw_pdb_dir.mkdir(parents=True, exist_ok=True)

    try:
        pdb_ids = list_pdb_ids_for_uniprot(accession)  # already best-resolution-first (server-side sort)
    except Exception as e:
        if show_progress:
            print(f"[templates] {accession}: RCSB lookup failed ({e}), skipping", flush=True)
        return []
    if not pdb_ids:
        if show_progress:
            print(f"[templates] {accession}: no RCSB structures", flush=True)
        return []

    kept: List[dict] = []
    for pdb_id in pdb_ids:
        try:
            meta = fetch_entry_metadata(pdb_id)
        except Exception:
            continue
        if meta.resolution is None or meta.resolution > resolution_cutoff:
            continue
        dest = raw_pdb_dir / f"{meta.pdb_id}.pdb"
        already_had_it = dest.exists()
        try:
            download_pdb(meta.pdb_id, dest)
        except Exception as e:
            if show_progress:
                print(f"[templates] {accession}: {meta.pdb_id} download failed ({e}), skipping", flush=True)
            continue
        kept.append({
            "pdb_id": meta.pdb_id, "resolution": meta.resolution, "method": meta.method,
            "title": meta.title, "pdb_path": str(dest),
        })
        if show_progress:
            status = "already downloaded" if already_had_it else "downloaded"
            print(
                f"[templates] {accession}: {meta.pdb_id} (resolution={meta.resolution}Å, "
                f"{meta.method}) {status} [{len(kept)} kept so far]", flush=True,
            )
    if show_progress and not kept:
        print(f"[templates] {accession}: {len(pdb_ids)} RCSB structure(s), none <= {resolution_cutoff}Å", flush=True)
    return kept
