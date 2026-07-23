"""Given a seed UniProt accession, find candidate similar proteins by Pfam/
InterPro family membership, ranked by global sequence identity to the seed.

Family selection: a protein is typically cross-referenced to several Pfam/
InterPro entries at different specificity levels -- e.g. human CDK20 carries
both `InterPro:IPR050108` ("CDK" family, 26 human reviewed members) and the
much broader `InterPro:IPR000719` ("Prot_kinase_dom", 480 members) and
`Pfam:PF00069` ("Pkinase", 344 members, effectively "has a kinase domain at
all"). Searching the broadest one would flood the candidate list with
proteins that share only a generic domain, not real similarity; the
narrowest ones can resolve to just the seed itself (1 member). This module
resolves the smallest cross-referenced family that still has more than one
member -- the most specific classification that actually contains other
proteins to compare against.

Known limitation (see README): this only finds proteins UniProt/InterPro
already classifies in the *same* family as the seed. A functionally related
but differently classified protein -- e.g. MAK relative to CDK20/CDK2, a
worked example in this project's own README -- will not surface here and
must be added to an explicit accession list instead.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional

from . import fetch
from .sequence import percent_identity


@dataclass
class FamilyChoice:
    database: str
    family_id: str
    entry_name: str
    member_count: int


def choose_family(entry: dict, *, taxon_id: Optional[int], show_progress: bool = True) -> Optional[FamilyChoice]:
    """The smallest Pfam/InterPro family cross-referenced on `entry` with
    more than one (reviewed, same-organism unless `taxon_id` is None)
    member. Ties broken by member count, then InterPro before Pfam (InterPro
    entries are generally the more curated/specific classification), then
    id. Returns None if every cross-referenced family resolves to just the
    seed itself (nothing to compare against)."""
    scored = []
    for db, fid, name in fetch.family_cross_references(entry):
        count = fetch.count_family_members(db, fid, taxon_id=taxon_id)
        scored.append((count, db, fid, name))
        if show_progress:
            print(f"[discover] candidate family {db}:{fid} ({name}): {count} member(s)", flush=True)
    eligible = [c for c in scored if c[0] > 1]
    if not eligible:
        return None
    count, db, fid, name = min(eligible, key=lambda c: (c[0], c[1] != "InterPro", c[2]))
    return FamilyChoice(database=db, family_id=fid, entry_name=name, member_count=count)


@dataclass
class Candidate:
    accession: str
    name: str
    length: int
    pct_identity: float


def discover(
    seed_accession: str, *, max_candidates: int = 20, any_organism: bool = False, show_progress: bool = True,
) -> dict:
    """Find and rank candidate similar proteins for `seed_accession`. Does
    not fetch or download anything beyond sequences/metadata -- this is a
    proposal step (writes `candidates.json`); the caller decides which
    candidates to actually pass to `pipeline.fetch_all`."""
    seed_accession = seed_accession.upper()
    seed_entry = fetch.fetch_uniprot_entry(seed_accession)
    seed_seq = fetch.fetch_uniprot_fasta(seed_accession)
    taxon_id = None if any_organism else fetch.organism_taxon_id(seed_entry)

    family = choose_family(seed_entry, taxon_id=taxon_id, show_progress=show_progress)
    if family is None:
        raise ValueError(
            f"{seed_accession}: every cross-referenced Pfam/InterPro family resolves to just this "
            f"protein itself -- nothing to compare against automatically. Pass an explicit accession "
            f"list to `dd_compare-fetch` instead."
        )
    if show_progress:
        print(
            f"[discover] {seed_accession}: using family {family.database}:{family.family_id} "
            f"({family.entry_name}), {family.member_count} member(s)", flush=True,
        )

    member_accessions = fetch.list_family_members(
        family.database, family.family_id, taxon_id=taxon_id, limit=max(family.member_count, max_candidates),
    )
    member_accessions = [a for a in member_accessions if a.upper() != seed_accession]

    candidates: List[Candidate] = []
    for i, acc in enumerate(member_accessions, start=1):
        try:
            seq = fetch.fetch_uniprot_fasta(acc)
            entry = fetch.fetch_uniprot_entry(acc)
        except Exception as e:
            if show_progress:
                print(f"[discover] ({i}/{len(member_accessions)}) {acc}: skipped ({e})", flush=True)
            continue
        pct = percent_identity(seed_seq, seq)
        name = fetch.protein_name(entry)
        candidates.append(Candidate(accession=acc, name=name, length=len(seq), pct_identity=pct))
        if show_progress:
            print(f"[discover] ({i}/{len(member_accessions)}) {acc} ({name}): {pct:.1f}% identity to seed", flush=True)

    candidates.sort(key=lambda c: c.pct_identity, reverse=True)
    candidates = candidates[:max_candidates]

    return {
        "seed_accession": seed_accession,
        "seed_name": fetch.protein_name(seed_entry),
        "family_database": family.database,
        "family_id": family.family_id,
        "family_entry_name": family.entry_name,
        "family_member_count": family.member_count,
        "candidates": [c.__dict__ for c in candidates],
    }
