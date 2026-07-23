"""Command-line entry points:
  dd_compare-fetch  ACC [ACC ...] -o out_dir
  dd_compare-fetch  --discover SEED_ACC -o out_dir   (writes candidates.json only)
  dd_compare-align  out_dir --reference ACC
  dd_compare-run    ACC [ACC ...] -o out_dir --reference ACC
"""
from __future__ import annotations

import argparse

from . import pipeline


def build_fetch_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_compare-fetch",
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
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")
    return parser


def _add_align_args(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--reference", default=None, help="Accession to superpose everything onto and detect the pocket on (default: the first accession fetched)")
    parser.add_argument("--pocket-rank", type=int, default=1, help="Druggability-ranked pocket to use on the reference (default: 1, top-ranked)")
    parser.add_argument(
        "--no-pdb-overlay", action="store_true",
        help="Skip looking up each protein's real RCSB structures. By default, when a protein has any, one "
             "(preferring a ligand-bound entry, else best resolution) is fetched and superposed onto the "
             "reference alongside its AlphaFold model, purely as an additional visual layer -- pocket detection "
             "and the cross-protein sequence/pocket mapping always stay anchored on the AlphaFold model.",
    )
    parser.add_argument(
        "--pdb-scan-cap", type=int, default=25,
        help="With PDB overlay enabled: how many resolution-ranked candidate structures to check for a bound "
             "ligand before falling back to the best-resolution one (default: 25) -- caps network/download cost "
             "for well-studied targets with hundreds of entries.",
    )
    parser.add_argument("--no-progress", action="store_true", help="Suppress the one-line-per-item progress output")


def build_align_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_compare-align",
        description="Detect the reference protein's druggable pocket, align every other protein's sequence to it "
                    "and map the pocket onto each, then superpose every protein's AlphaFold model onto the reference.",
    )
    parser.add_argument("out_dir", help="Directory previously populated by dd_compare-fetch")
    _add_align_args(parser)
    return parser


def build_run_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dd_compare-run",
        description="dd_compare-fetch followed by dd_compare-align in one step (explicit accession list only -- "
                    "run dd_compare-fetch --discover separately first if you want candidate suggestions).",
    )
    parser.add_argument("accessions", nargs="+", help="UniProt accessions, e.g. Q8IZL9 P24941 P20794")
    parser.add_argument("-o", "--out-dir", required=True, help="Output directory")
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
        raise SystemExit("dd_compare-fetch: pass at least 2 accessions to compare, or --discover SEED_ACC")
    manifest = pipeline.fetch_all(args.accessions, args.out_dir, show_progress=not args.no_progress)
    print(f"\n[done] {len(manifest['proteins'])} protein(s) -> {args.out_dir}")


def main_align(argv=None) -> None:
    args = build_align_parser().parse_args(argv)
    report = pipeline.analyze(
        args.out_dir, reference=args.reference, pocket_rank=args.pocket_rank, show_progress=not args.no_progress,
        pdb_overlay=not args.no_pdb_overlay, pdb_scan_cap=args.pdb_scan_cap,
    )
    _print_report(report, args.out_dir)


def main_run(argv=None) -> None:
    args = build_run_parser().parse_args(argv)
    if len(args.accessions) < 2:
        raise SystemExit("dd_compare-run: pass at least 2 accessions to compare")
    pipeline.fetch_all(args.accessions, args.out_dir, show_progress=not args.no_progress)
    report = pipeline.analyze(
        args.out_dir, reference=args.reference, pocket_rank=args.pocket_rank, show_progress=not args.no_progress,
        pdb_overlay=not args.no_pdb_overlay, pdb_scan_cap=args.pdb_scan_cap,
    )
    _print_report(report, args.out_dir)


def _print_report(report: dict, out_dir: str) -> None:
    print(f"\n[reference: {report['reference']}] pocket: {len(report['pocket']['residues'])} residue(s)")
    for p in report["proteins"]:
        rmsd = f"{p['rmsd']:.3f} ({p['n_aligned_atoms']} atoms)" if p["rmsd"] is not None else f"SKIPPED: {p['align_error']}"
        n_noncons = sum(1 for c in p["pocket_comparison"] if c["conservation"] == "non-conservative")
        print(f"  {p['accession']:<10} {p['name']:<40.40} identity={p['pct_identity']:>5.1f}%  rmsd={rmsd}  non-conservative-pocket-residues={n_noncons}")
    print(f"\n[done] report -> {out_dir}/report.json")
