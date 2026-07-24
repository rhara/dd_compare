"""Command-line entry points:
  dd_idea-fetch  ACC [ACC ...] -o out_dir
  dd_idea-fetch  --discover SEED_ACC -o out_dir   (writes candidates.json only)
  dd_idea-align  out_dir --reference ACC
  dd_idea-run    ACC [ACC ...] -o out_dir --reference ACC
  dd_idea-search QUERY -o out_dir   (UniProt accession, ChEMBL target ID, or a raw sequence -- table only)
  dd_idea-search --fetch ACC [ACC ...] -o out_dir   (or --fetch-all -- downloads templates for hits.json rows)
"""
from __future__ import annotations

import argparse
import csv

from . import pipeline, search as search_module


def _add_pdb_fetch_args(parser: argparse.ArgumentParser) -> None:
    """Real-RCSB-structure lookup/selection flags -- genuinely fetch-time
    network work (see `pipeline.fetch_all`), so shared by `dd_idea-fetch`
    and `dd_idea-run` only, not `dd_idea-align` (which just reuses
    whatever was already fetched, with no network access -- re-running
    align never re-hits RCSB, so it has nothing to configure here)."""
    parser.add_argument(
        "--no-pdb-overlay", action="store_true",
        help="Skip looking up each protein's real RCSB structures. By default, when a protein has any, up to "
             "--pdb-max-structures of them (preferring distinct ligand-bound entries, else best resolution) are "
             "fetched for later superposition onto the reference alongside its AlphaFold model, purely as an "
             "additional visual layer -- pocket detection and the cross-protein sequence/pocket mapping always "
             "stay anchored on the AlphaFold model.",
    )
    parser.add_argument(
        "--pdb-max-structures", type=int, default=3,
        help="With PDB overlay enabled: how many distinct-ligand real structures to fetch per protein (default: 3). "
             "Re-running fetch with a larger value cheaply extends what's already cached on disk.",
    )
    parser.add_argument(
        "--pdb-scan-cap", type=int, default=25,
        help="With PDB overlay enabled: how many resolution-ranked candidate structures to check for a bound "
             "ligand, at most, before giving up on finding --pdb-max-structures of them and falling back to the "
             "best-resolution one (default: 25) -- caps network/download cost for well-studied targets with "
             "hundreds of entries.",
    )
    parser.add_argument(
        "--pdb-resolution-cutoff", type=float, default=2.0,
        help="With PDB overlay enabled: worst resolution (Angstrom, lower is better) a real structure may have to "
             "be considered at all (default: 2.0) -- entries worse than this, or with no reported resolution "
             "(e.g. NMR), are skipped before download, for both the ligand-bound and best-resolution-fallback "
             "paths.",
    )


def build_fetch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_idea-fetch",
        description="Download the canonical sequence + AlphaFold DB model for two or more UniProt accessions, "
                    "or (with --discover) propose candidate similar proteins for one seed accession.",
    )
    parser.add_argument("accessions", nargs="*", help="UniProt accessions, e.g. Q8IZL9 P24941 P20794")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument(
        "--discover", metavar="SEED_ACC", default=None,
        help="Instead of fetching, propose candidate similar proteins for this seed accession "
             "(Pfam/InterPro family + sequence-identity ranking) and write candidates.json. "
             "Review the ranked list, then pass the accessions you want as this command's normal arguments.",
    )
    parser.add_argument("--max-candidates", type=int, default=20, help="With --discover: candidates to keep (default: 20)")
    parser.add_argument(
        "--any-organism", action="store_true",
        help="With --discover: don't restrict candidates to the seed's own organism",
    )
    _add_pdb_fetch_args(parser)
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")
    return parser


def _add_align_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reference", default=None, help="Accession to superpose everything onto and detect the pocket on (default: the first accession fetched)")
    parser.add_argument("--pocket-rank", type=int, default=1, help="Druggability-ranked pocket to use on the reference (default: 1, top-ranked)")
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")


def build_align_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_idea-align",
        description="Detect the reference protein's druggable pocket, align every other protein's sequence to it "
                    "and map the pocket onto each, then superpose every protein's AlphaFold model (and any real "
                    "RCSB structures dd_idea-fetch already selected) onto the reference.",
    )
    parser.add_argument("out_dir", help="Directory previously populated by dd_idea-fetch")
    _add_align_args(parser)
    return parser


def build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_idea-run",
        description="dd_idea-fetch followed by dd_idea-align in one step (explicit accession list only -- "
                    "run dd_idea-fetch --discover separately first if you want candidate suggestions).",
    )
    parser.add_argument("accessions", nargs="+", help="UniProt accessions, e.g. Q8IZL9 P24941 P20794")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    _add_pdb_fetch_args(parser)
    _add_align_args(parser)
    return parser


def main_fetch(argv=None) -> None:
    args = build_fetch_parser().parse_args(argv)
    if args.discover:
        result = pipeline.discover(
            args.discover, args.out_dir, max_candidates=args.max_candidates,
            any_organism=args.any_organism, show_progress=not args.no_progress,
        )
        print(
            f"\n[done] {len(result['candidates'])} candidate(s) (family {result['family_database']}:"
            f"{result['family_id']} {result['family_entry_name']}) -> {args.out_dir}/candidates.json"
        )
        for c in result["candidates"]:
            print(f"  {c['accession']:<10} {c['pct_identity']:>5.1f}%  {c['name']}")
        return

    if len(args.accessions) < 2:
        raise SystemExit("dd_idea-fetch: pass at least 2 accessions to compare, or --discover SEED_ACC")
    manifest = pipeline.fetch_all(
        args.accessions, args.out_dir, show_progress=not args.no_progress,
        pdb_overlay=not args.no_pdb_overlay, pdb_scan_cap=args.pdb_scan_cap,
        pdb_max_structures=args.pdb_max_structures, pdb_resolution_cutoff=args.pdb_resolution_cutoff,
    )
    print(f"\n[done] {len(manifest['proteins'])} protein(s) -> {args.out_dir}")


def main_align(argv=None) -> None:
    args = build_align_parser().parse_args(argv)
    report = pipeline.analyze(
        args.out_dir, reference=args.reference, pocket_rank=args.pocket_rank, show_progress=not args.no_progress,
    )
    _print_report(report, args.out_dir)


def main_run(argv=None) -> None:
    args = build_run_parser().parse_args(argv)
    if len(args.accessions) < 2:
        raise SystemExit("dd_idea-run: pass at least 2 accessions to compare")
    pipeline.fetch_all(
        args.accessions, args.out_dir, show_progress=not args.no_progress,
        pdb_overlay=not args.no_pdb_overlay, pdb_scan_cap=args.pdb_scan_cap,
        pdb_max_structures=args.pdb_max_structures, pdb_resolution_cutoff=args.pdb_resolution_cutoff,
    )
    report = pipeline.analyze(
        args.out_dir, reference=args.reference, pocket_rank=args.pocket_rank, show_progress=not args.no_progress,
    )
    _print_report(report, args.out_dir)


def _print_report(report: dict, out_dir: str) -> None:
    print(f"\n[reference: {report['reference']}] pocket: {len(report['pocket']['residues'])} residue(s)")
    for p in report["proteins"]:
        rmsd = f"{p['rmsd']:.3f} ({p['n_aligned_atoms']} atoms)" if p["rmsd"] is not None else f"SKIPPED: {p['align_error']}"
        n_noncons = sum(1 for c in p["pocket_comparison"] if c["conservation"] == "non-conservative")
        print(f"  {p['accession']:<10} {p['name']:<40.40} identity={p['pct_identity']:>5.1f}%  rmsd={rmsd}  non-conservative-pocket-residues={n_noncons}")
    print(f"\n[done] report -> {out_dir}/report.json")


def build_search_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_idea-search",
        description="Resolve QUERY (a UniProt accession, a ChEMBL target ID, or a raw amino-acid sequence) to a "
                    "protein and BLASTP it against Swiss-Prot to find real sequence-similar proteins (independent "
                    "of Pfam/InterPro family classification -- see --discover for that), printing a reviewable "
                    "table (family/gene/organism/length/identity) -- no downloads yet. Once you've looked at the "
                    "table, --fetch/--fetch-all downloads AlphaFold models + RCSB structures for just the "
                    "accessions worth it.",
    )
    parser.add_argument("query", nargs="?", default=None, help="UniProt accession (e.g. Q8IZL9), ChEMBL target ID (e.g. CHEMBL301), or a raw sequence -- required unless --fetch/--fetch-all")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
    parser.add_argument("--evalue", type=float, default=1e-10, help="BLASTP e-value threshold (default: 1e-10)")
    parser.add_argument("--any-organism", action="store_true", help="Search all organisms instead of the default Homo sapiens-only restriction")
    parser.add_argument(
        "--max-hits", type=int, default=100,
        help="Keep at most this many top BLAST hits, ranked by %%identity (default: 100, matching NCBI's own "
             "hitlist_size so nothing already returned is silently dropped)",
    )
    parser.add_argument(
        "--fetch", metavar="ACC", nargs="+", default=None,
        help="Instead of a new search, download AlphaFold model (if it's the seed) + RCSB structures for these "
             "accessions from an existing hits.json in --out-dir (run a plain dd_idea-search QUERY -o DIR first).",
    )
    parser.add_argument("--fetch-all", action="store_true", help="Like --fetch, but for every row in hits.json")
    parser.add_argument(
        "--resolution-cutoff", type=float, default=2.0,
        help="With --fetch/--fetch-all: worst resolution (Angstrom, lower is better) a real structure may have to "
             "be kept (default: 2.0). Every structure meeting this bar is kept -- no cap, no ligand preference "
             "(unlike dd_idea-fetch's --pdb-resolution-cutoff).",
    )
    parser.add_argument("--summary-format", choices=["table", "csv", "markdown"], default="table", help="table (stdout only, default) also written as csv/markdown to out_dir")
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")
    return parser


def _summary_headers_row(r: dict) -> list:
    evalue_str = "-" if r["role"] == "seed" else f"{r['evalue']:.2g}"
    n_templates = "-" if r["pdb_structures"] is None else str(len(r["pdb_structures"]))
    return [
        str(r["accession"] or "(sequence)"), r["gene"], r["organism"], r["family"], str(r["length"]),
        f"{r['pct_identity']:.1f}", evalue_str, n_templates,
    ]


def _print_search_table(rows: list) -> None:
    headers = ["Accession", "Gene", "Organism", "Family", "Length", "%Id", "E-value", "#Templates"]
    table_rows = [_summary_headers_row(r) for r in rows]
    widths = [max(len(headers[i]), *(len(row[i]) for row in table_rows)) for i in range(len(headers))]
    widths[3] = min(widths[3], 64)  # cap Family so one BLOSUM-style long comment doesn't blow out the table
    fmt = "  ".join(f"{{:<{w}}}" for w in widths)
    print(fmt.format(*headers))
    print(fmt.format(*["-" * w for w in widths]))
    for row in table_rows:
        row = list(row)
        row[3] = row[3][:64]
        print(fmt.format(*row))


def _write_search_summary(rows: list, out_dir: str, fmt: str) -> None:
    headers = ["Accession", "Gene", "Organism", "Family", "Length", "%Id", "E-value", "#Templates"]
    table_rows = [_summary_headers_row(r) for r in rows]
    if fmt == "csv":
        dest = f"{out_dir}/hits_summary.csv"
        with open(dest, "w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(headers)
            writer.writerows(table_rows)
    else:
        dest = f"{out_dir}/hits_summary.md"
        lines = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
        lines += ["| " + " | ".join(row) + " |" for row in table_rows]
        with open(dest, "w") as fh:
            fh.write("\n".join(lines) + "\n")
    print(f"\n-> {dest}")


def main_search(argv=None) -> None:
    args = build_search_parser().parse_args(argv)
    show_progress = not args.no_progress

    if args.fetch_all or args.fetch:
        if args.query:
            print(f"[note] ignoring QUERY ({args.query!r}) -- --fetch/--fetch-all operates on the existing hits.json in {args.out_dir}")
        accessions = "all" if args.fetch_all else args.fetch
        result = search_module.fetch_templates(
            args.out_dir, accessions, resolution_cutoff=args.resolution_cutoff, show_progress=show_progress,
        )
        rows = result["hits"]
        print()
        _print_search_table(rows)
        if args.summary_format in ("csv", "markdown"):
            _write_search_summary(rows, args.out_dir, args.summary_format)
        n_structures = sum(len(r["pdb_structures"]) for r in rows if r["pdb_structures"] is not None)
        print(f"\n[done] templates fetched -> {n_structures} structure(s) total -> {args.out_dir}/hits.json")
        return

    if not args.query:
        raise SystemExit("dd_idea-search: QUERY is required unless --fetch/--fetch-all is given")
    result = search_module.search(
        args.query, args.out_dir, evalue=args.evalue, any_organism=args.any_organism,
        max_hits=args.max_hits, show_progress=show_progress,
    )
    rows = result["hits"]
    print()
    _print_search_table(rows)
    if args.summary_format in ("csv", "markdown"):
        _write_search_summary(rows, args.out_dir, args.summary_format)
    print(
        f"\n[done] {len(rows)} protein(s) (seed + hits) -> {args.out_dir}/hits.json -- "
        f"no downloads yet, review the table then run with --fetch ACC [ACC ...] or --fetch-all"
    )
