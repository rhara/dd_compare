"""`dd_idea-search`: the general entry point into `dd_idea` -- given a
UniProt accession, a ChEMBL target ID, or a raw amino-acid sequence
pasted directly, resolve it to a protein, BLAST it against Swiss-Prot to
find real sequence-similar proteins (see `blast/`'s docstring for why
this exists alongside `similarity.py`'s Pfam/InterPro-based `--discover`),
and build a reviewable table (UniProt family, gene, organism, length,
%identity, e-value) -- deliberately *without* downloading anything yet.

Atomic by design: `search()` only ever talks to UniProt + NCBI BLAST
(small, fast requests) and writes `hits.json`. Two further, separate,
explicit steps enrich specific rows of that same `hits.json` only for
accessions the caller actually asks for, after looking at the table
`search()` produced -- `fetch_templates()` downloads AlphaFold models and
RCSB structures (a protein with hundreds of RCSB entries, e.g. CDK2's 512,
should never be downloaded by accident just because it BLAST-matched the
seed), and `annotate_chembl_activity()` counts each accession's ChEMBL
bioactivity records (cheap -- counts only, no data itself) to gauge how
much SAR data exists for it.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Union

from . import alphafold, blast, chembl, rcsb, uniprot

_CHEMBL_RE = re.compile(r"^CHEMBL\d+$", re.IGNORECASE)
# UniProt canonical accession: 6 chars (e.g. Q8IZL9, P24941) or 10 chars.
_UNIPROT_ACC_RE = re.compile(
    r"^([OPQ][0-9][A-Z0-9]{3}[0-9])$|^([A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$"
)
_AA_LETTERS = set("ACDEFGHIKLMNPQRSTVWYXBZUO")


def classify_query(query: str) -> str:
    """`"chembl"`, `"uniprot"`, or `"sequence"` -- pure string
    classification, no network access. Raises `ValueError` (listing the
    three accepted forms) for anything else, e.g. too-short garbage that's
    neither a recognizable accession/id nor a plausible sequence."""
    q = query.strip().upper()
    if _CHEMBL_RE.match(q):
        return "chembl"
    if _UNIPROT_ACC_RE.match(q):
        return "uniprot"
    if len(q) >= 15 and all(c in _AA_LETTERS for c in q):
        return "sequence"
    raise ValueError(
        f"{query!r} isn't a recognizable UniProt accession (e.g. Q8IZL9), ChEMBL target ID (e.g. CHEMBL301), "
        f"or amino-acid sequence (>=15 residues, one-letter codes)."
    )


@dataclass
class SeedInfo:
    query_kind: str
    accession: Optional[str]
    canonical_seq: str
    entry: Optional[dict]
    fasta_path: Optional[str]


def resolve_seed(query: str, out_dir: Union[str, Path], *, show_progress: bool = True) -> SeedInfo:
    """Resolve `query` to a canonical sequence (+ UniProt entry, when an
    accession is known) -- the minimum needed to BLAST and to fill in the
    table's metadata columns. Does **not** fetch an AlphaFold model; that
    only happens in `fetch_templates`, once a specific accession is
    actually wanted. Branches on `classify_query`: a ChEMBL target ID is
    resolved to a UniProt accession first (`chembl.resolve_chembl_target_to_uniprot`),
    then handled exactly like a direct UniProt-accession query. A raw
    sequence is used as-is, with `accession=None`/`entry=None` -- there's
    no UniProt entry to fetch metadata from for an accession-less
    sequence."""
    out_dir = Path(out_dir)
    kind = classify_query(query)

    if kind == "sequence":
        seq = query.strip().upper()
        if show_progress:
            print(f"[seed] raw sequence ({len(seq)} aa), no UniProt accession -- skipping metadata fetch", flush=True)
        return SeedInfo(query_kind=kind, accession=None, canonical_seq=seq, entry=None, fasta_path=None)

    accession = query.strip().upper()
    if kind == "chembl":
        accession = chembl.resolve_chembl_target_to_uniprot(accession)
        if show_progress:
            print(f"[seed] {query.strip().upper()}: resolved to UniProt {accession}", flush=True)

    out_dir.mkdir(parents=True, exist_ok=True)
    entry = uniprot.fetch_uniprot_entry(accession)

    fasta_dest = out_dir / f"{accession}.fasta"
    if fasta_dest.exists():
        canonical = "".join(fasta_dest.read_text().splitlines()[1:])
        if show_progress:
            print(f"[seed] {accession}: canonical sequence already downloaded, skipping", flush=True)
    else:
        canonical = uniprot.fetch_uniprot_fasta(accession)
        fasta_dest.write_text(f">{accession}\n{canonical}\n")
        if show_progress:
            print(f"[seed] {accession}: canonical sequence ({len(canonical)} aa)", flush=True)

    return SeedInfo(query_kind=kind, accession=accession, canonical_seq=canonical, entry=entry, fasta_path=str(fasta_dest))


def _metadata_row(accession: Optional[str], entry: Optional[dict], *, pct_identity: float, evalue: Optional[float], length: Optional[int] = None, role: str) -> dict:
    """A `hits.json` row with UniProt metadata filled in but
    `pdb_structures: None` (not yet fetched via `fetch_templates`) and
    `chembl_targets: None` (not yet annotated via `annotate_chembl_activity`)
    -- the shape `search()` writes for every row, before either later,
    optional enrichment step has touched anything."""
    if entry is None:
        return {
            "accession": accession, "gene": "-", "organism": "-", "family": "(no UniProt accession)",
            "length": length, "pct_identity": pct_identity, "evalue": evalue,
            "afdb_path": None, "pdb_structures": None, "chembl_targets": None, "role": role,
        }
    return {
        "accession": accession, "gene": uniprot.gene_name(entry), "organism": uniprot.organism_name(entry),
        "family": uniprot.family_string(entry), "length": uniprot.canonical_length(entry),
        "pct_identity": pct_identity, "evalue": evalue, "afdb_path": None, "pdb_structures": None,
        "chembl_targets": None, "role": role,
    }


def search(
    query: str, out_dir: Union[str, Path], *,
    evalue: float = 1e-10, any_organism: bool = False, max_hits: int = 100, show_progress: bool = True,
) -> dict:
    """Resolve `query` to a seed protein, BLAST it against Swiss-Prot, and
    fetch UniProt metadata for the seed (if it has an accession) plus
    every hit. Writes and returns `{out_dir}/hits.json` -- every row's
    `pdb_structures` is `None` (distinct from `[]`, "fetched, found
    none"). No AlphaFold/RCSB downloads happen here; call
    `fetch_templates` next for the accessions worth pulling structures
    for."""
    out_dir = Path(out_dir)
    seed = resolve_seed(query, out_dir, show_progress=show_progress)

    xml_text = blast.run_blastp(
        seed.canonical_seq, out_dir, evalue=evalue, any_organism=any_organism, show_progress=show_progress,
    )
    hits = blast.parse_blast_hits(xml_text, exclude_accession=seed.accession, max_hits=max_hits)
    if show_progress:
        print(f"[blast] {len(hits)} Swiss-Prot hit(s) kept", flush=True)

    rows: List[dict] = [_metadata_row(
        seed.accession, seed.entry, pct_identity=100.0, evalue=None, length=len(seed.canonical_seq), role="seed",
    )]
    for i, hit in enumerate(hits, start=1):
        if show_progress:
            print(f"[uniprot] ({i}/{len(hits)}) {hit.accession}: fetching metadata...", flush=True)
        entry = uniprot.fetch_uniprot_entry(hit.accession)
        rows.append(_metadata_row(hit.accession, entry, pct_identity=hit.pct_identity, evalue=hit.evalue, role="blast_hit"))

    result = {
        "query": query, "query_kind": seed.query_kind, "seed_accession": seed.accession,
        "blast": {
            "database": "swissprot", "evalue_threshold": evalue,
            "organism_restricted": not any_organism, "max_hits": max_hits,
        },
        "hits": rows,
    }
    (out_dir / "hits.json").write_text(json.dumps(result, indent=2))
    return result


def _load_hits_json(out_dir: Path) -> dict:
    hits_path = out_dir / "hits.json"
    if not hits_path.exists():
        raise FileNotFoundError(
            f"{hits_path} not found -- run `dd_idea-search QUERY -o {out_dir}` first to build the table."
        )
    return json.loads(hits_path.read_text())


def _resolve_targets(result: dict, accessions: Union[List[str], str]) -> "dict[str, dict]":
    """`{accession: row}` for the requested `accessions` (or every row, for
    the literal string `"all"`) out of an already-loaded `hits.json`.
    Raises `ValueError` naming the offending accession(s) if any aren't
    present -- shared by `fetch_templates` and `annotate_chembl_activity`,
    the two `hits.json`-enrichment steps."""
    rows_by_acc = {r["accession"]: r for r in result["hits"] if r["accession"] is not None}
    if accessions == "all":
        return rows_by_acc
    unknown = [a for a in accessions if a.upper() not in rows_by_acc]
    if unknown:
        raise ValueError(f"unknown accession(s) {unknown} -- hits.json has: {sorted(rows_by_acc)}")
    return {a.upper(): rows_by_acc[a.upper()] for a in accessions}


def _gene_subdir(acc: str, row: dict) -> str:
    """The `raw_pdb/` subdirectory a row's structures download into --
    its gene symbol (e.g. `CDK2`), falling back to the accession itself
    when no gene name is on record (e.g. a raw-sequence seed with no
    UniProt entry at all)."""
    gene = row.get("gene")
    return gene if gene and gene != "-" else acc


def fetch_templates(
    out_dir: Union[str, Path], accessions: Union[List[str], str], *,
    resolution_cutoff: float = 2.0, show_progress: bool = True,
) -> dict:
    """Download AlphaFold model (seed row only) + every RCSB structure at
    or better than `resolution_cutoff` for the given `accessions` (or the
    literal string `"all"` for every row) out of an existing
    `{out_dir}/hits.json` (written by a prior `search()` call). Updates
    those rows in place and rewrites `hits.json`. Each accession's RCSB
    structures land under their own `raw_pdb/{gene}/` subdirectory (see
    `_gene_subdir`), not flat in `raw_pdb/` -- a multi-accession fetch can
    mean hundreds of files from unrelated proteins otherwise."""
    out_dir = Path(out_dir)
    result = _load_hits_json(out_dir)
    targets = _resolve_targets(result, accessions)

    for i, (acc, row) in enumerate(targets.items(), start=1):
        if show_progress:
            print(f"[fetch] ({i}/{len(targets)}) {acc}: fetching templates...", flush=True)
        if row["role"] == "seed":
            afdb_dest = out_dir / "raw" / f"{acc}_AFDB.pdb"
            try:
                alphafold.download_afdb(acc, afdb_dest)
                row["afdb_path"] = str(afdb_dest)
                if show_progress:
                    print(f"[fetch] {acc}: AlphaFold DB model -> {afdb_dest.name}", flush=True)
            except Exception as e:
                if show_progress:
                    print(f"[fetch] {acc}: AlphaFold DB model unavailable ({e})", flush=True)
        row["pdb_structures"] = rcsb.list_all_structures_at_resolution(
            acc, out_dir, resolution_cutoff=resolution_cutoff, subdir=_gene_subdir(acc, row), show_progress=show_progress,
        )

    (out_dir / "hits.json").write_text(json.dumps(result, indent=2))
    return result


def annotate_chembl_activity(
    out_dir: Union[str, Path], accessions: Union[List[str], str], *,
    assay_types: tuple = chembl.DEFAULT_ASSAY_TYPES, show_progress: bool = True,
) -> dict:
    """For the given `accessions` (or `"all"`) out of an existing
    `hits.json`, resolve every ChEMBL `SINGLE PROTEIN` target and count its
    bioactivity records (`chembl.count_activities` -- binding assays with a
    pChEMBL value by default, same filter `dd_chembl` itself uses for QSAR
    training data). Sets each row's `chembl_targets` to a list of
    `{target_chembl_id, pref_name, organism, n_activities}` (`[]` if
    ChEMBL has no single-protein target for that accession at all -- not
    the same as `None`, "not yet checked"). Cheap: no bioactivity data
    itself is downloaded, only counts."""
    out_dir = Path(out_dir)
    result = _load_hits_json(out_dir)
    targets = _resolve_targets(result, accessions)

    for i, (acc, row) in enumerate(targets.items(), start=1):
        if show_progress:
            print(f"[chembl] ({i}/{len(targets)}) {acc}: resolving ChEMBL target(s)...", flush=True)
        try:
            chembl_targets = chembl.resolve_uniprot_to_chembl_targets(acc)
        except Exception as e:
            if show_progress:
                print(f"[chembl] {acc}: target resolution failed ({e}), skipping", flush=True)
            continue
        if not chembl_targets:
            if show_progress:
                print(f"[chembl] {acc}: no ChEMBL single-protein target", flush=True)
            row["chembl_targets"] = []
            continue
        annotated = []
        for t in chembl_targets:
            try:
                n = chembl.count_activities(t["target_chembl_id"], assay_types=assay_types)
            except Exception as e:
                if show_progress:
                    print(f"[chembl] {acc}: {t['target_chembl_id']} activity count failed ({e})", flush=True)
                continue
            annotated.append({**t, "n_activities": n})
            if show_progress:
                print(f"[chembl] {acc}: {t['target_chembl_id']} ({t['pref_name']}): {n} activities", flush=True)
        row["chembl_targets"] = annotated

    (out_dir / "hits.json").write_text(json.dumps(result, indent=2))
    return result
