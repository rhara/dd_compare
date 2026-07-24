"""NCBI QBLAST submission (the network/caching half of `blast/`; see
`blast/parse.py` for turning a result into `BlastHit`s -- kept separate so
re-parsing an already-cached XML never needs this module at all).

Talks to NCBI's BLAST URL API directly (rather than through
`Bio.Blast.NCBIWWW.qblast`) so the RID (Request ID, NCBI's acknowledgement
that it received and queued the search) and each poll's status can be
printed as they happen -- `qblast` receives both internally but never
surfaces them, so a slow or stuck search is otherwise silent until
Biopython's own 10-minutes-and-counting warning fires."""
from __future__ import annotations

import re
import socket
import time
from contextlib import contextmanager
from pathlib import Path
from typing import Union
from urllib.parse import urlencode
from urllib.request import Request, urlopen

_NCBI_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"
_USER_AGENT = "dd_idea (github.com/rhara/dd_idea)"


@contextmanager
def _ipv4_only():
    """blast.ncbi.nlm.nih.gov resolves to both an IPv6 and an IPv4 address;
    on networks where IPv6 is routed but silently black-holed, `urllib`
    (unlike curl's Happy Eyeballs) tries the IPv6 address first and blocks
    for the OS's full TCP connect timeout -- tens of seconds -- before ever
    trying IPv4. Restrict resolution to IPv4 for the duration of the call."""
    original = socket.getaddrinfo

    def _getaddrinfo_ipv4(host, port, family=0, type=0, proto=0, flags=0):
        return original(host, port, socket.AF_INET, type, proto, flags)

    socket.getaddrinfo = _getaddrinfo_ipv4
    try:
        yield
    finally:
        socket.getaddrinfo = original


def _post(params: dict) -> str:
    request = Request(_NCBI_URL, urlencode(params).encode(), {"User-Agent": _USER_AGENT})
    with _ipv4_only(), urlopen(request, timeout=60) as handle:
        return handle.read().decode()


def run_blastp(
    sequence: str, out_dir: Union[str, Path], *,
    evalue: float = 1e-10, any_organism: bool = False, hitlist_size: int = 100, show_progress: bool = True,
) -> str:
    """BLASTP `sequence` against Swiss-Prot via NCBI's BLAST URL API. Cached:
    skipped (just re-read) if `{out_dir}/raw_blast/blastp_swissprot.xml`
    already exists -- a single search already takes NCBI several minutes."""
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
        print(f"[blast] submitting blastp vs. swissprot ({scope}, evalue<{evalue}) to NCBI...", flush=True)

    put_params = {
        "CMD": "Put", "PROGRAM": "blastp", "DATABASE": "swissprot",
        "QUERY": sequence, "EXPECT": evalue, "HITLIST_SIZE": hitlist_size,
    }
    if entrez_query:
        put_params["ENTREZ_QUERY"] = entrez_query
    put_response = _post(put_params)

    rid_match = re.search(r"RID = (\S+)", put_response)
    if not rid_match:
        raise RuntimeError(
            "NCBI did not return a request ID (RID) -- the submission itself failed; "
            "check network connectivity and try again"
        )
    rid = rid_match.group(1)
    rtoe_match = re.search(r"RTOE = (\d+)", put_response)
    rtoe = int(rtoe_match.group(1)) if rtoe_match else None

    if show_progress:
        eta = f", estimated ~{rtoe}s" if rtoe else ""
        print(
            f"[blast] NCBI acknowledged the request -- RID={rid}{eta}. Check status anytime at "
            f"https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID={rid}",
            flush=True,
        )

    # NCBI usage guidelines: don't poll more than once/minute per RID. Use
    # RTOE as the first wait, but cap it at that same 60s -- RTOE is a rough
    # server-load-dependent estimate and can come back in the thousands of
    # seconds on a busy queue, which would otherwise leave the first status
    # check (and thus the first sign of life after submission) that far away.
    poll_params = {"CMD": "Get", "FORMAT_OBJECT": "SearchInfo", "RID": rid}
    delay = min(max(rtoe or 20, 20), 60)
    start = time.time()
    while True:
        time.sleep(delay)
        delay = 60
        status_response = _post(poll_params)
        status_match = re.search(r"Status=(\S+)", status_response)
        status = status_match.group(1) if status_match else "UNKNOWN"
        elapsed = int(time.time() - start)
        if show_progress:
            print(f"[blast] ({elapsed}s elapsed) NCBI status: {status}", flush=True)
        if status == "READY":
            break
        if status == "FAILED":
            raise RuntimeError(f"NCBI reported that BLAST search {rid} failed")
        if status == "UNKNOWN":
            raise RuntimeError(f"NCBI has no record of RID {rid} (expired, or the RID itself was invalid)")

    xml_text = _post({"CMD": "Get", "FORMAT_TYPE": "XML", "RID": rid})
    dest.write_text(xml_text)
    if show_progress:
        print(f"[blast] done -> {dest}", flush=True)
    return xml_text
