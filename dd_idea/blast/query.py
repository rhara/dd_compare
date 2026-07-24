"""NCBI QBLAST submission (the network/caching half of `blast/`; see
`blast/parse.py` for turning a result into `BlastHit`s -- kept separate so
re-parsing an already-cached XML never needs this module at all)."""
from __future__ import annotations

from pathlib import Path
from typing import Union

from Bio.Blast import NCBIWWW


def run_blastp(
    sequence: str, out_dir: Union[str, Path], *,
    evalue: float = 1e-10, any_organism: bool = False, hitlist_size: int = 100, show_progress: bool = True,
) -> str:
    """BLASTP `sequence` against Swiss-Prot via NCBI QBLAST. Cached:
    skipped (just re-read) if `{out_dir}/raw_blast/blastp_swissprot.xml`
    already exists -- a single search already takes NCBI several
    minutes."""
    raw_dir = Path(out_dir) / "raw_blast"
    raw_dir.mkdir(parents=True, exist_ok=True)
    dest = raw_dir / "blastp_swissprot.xml"
    if dest.exists():
        if show_progress:
            print(f"[blast] blastp vs. swissprot already run, reusing {dest}", flush=True)
        return dest.read_text()

    entrez_query = None if any_organism else "Homo sapiens[Organism]"
    if show_progress:
        scope = "any organism" if any_organism else "Homo sapiens only"
        print(
            f"[blast] submitting blastp vs. swissprot ({scope}, evalue<{evalue}) to NCBI -- "
            f"this typically takes a few minutes...", flush=True,
        )
    handle = NCBIWWW.qblast(
        "blastp", "swissprot", sequence, expect=evalue, entrez_query=entrez_query, hitlist_size=hitlist_size,
    )
    xml_text = handle.read()
    handle.close()
    dest.write_text(xml_text)
    if show_progress:
        print(f"[blast] done -> {dest}", flush=True)
    return xml_text
