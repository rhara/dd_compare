"""Ranking `hits.json` rows by combined usefulness -- sequence identity to
the seed, family relatedness, RCSB template availability, and ChEMBL
bioactivity coverage. Pure computation over data `search.py`'s other
functions already gathered (no network access, no new fetching): a
`dd_idea-search --rank` pass is instant no matter how much `--fetch`/
`--chembl-activity` work preceded it.

Each of the four signals is converted to an ordinal class before
combining -- *not* multiplied as raw numbers -- because their raw scales
are wildly different and skewed in different ways (e.g. ChEMBL activity
counts span 0-7000+ while identity spans a compact 25-45% band; template
counts are zero for roughly half of any typical hit set). Multiplying raw
magnitudes would let whichever metric happens to have the largest numbers
dominate the ranking regardless of what it actually means biologically.
Identity and family use 1-5 classes -- family in particular is a genuinely
discrete signal (exact subfamily match / one level short / superfamily
only / no match), so finer bins would just split ties, not add real
information. Templates and activity use a wider 1-20 by default -- their
raw ranges are wide enough (hundreds of templates, thousands of
activities) that 5 quantile bins collapsed real differences into ties
(e.g. CDK2's 522 templates/3015 activities and CDK9's 28/2051 both
landing in template/activity class 5 despite being nowhere close). The
composite score is the product of the four classes (range 1-10000 with
the defaults): a hit weak on any single axis is penalized
multiplicatively, not just averaged away by the other three.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Sequence, Union

from .search import _load_hits_json

N_CLASSES = 5
N_COUNT_CLASSES = 20  # templates_class/activity_class default -- see module docstring


def classify_zero_inflated(value: float, all_values: Sequence[float], *, n_classes: int = N_CLASSES) -> int:
    """Class 1 for `value == 0` always (a whole class of its own, not
    folded into the bottom quantile -- appropriate when a large fraction
    of `all_values` are exactly zero, e.g. RCSB template counts: about
    half of a typical BLAST hit set has none at all). Nonzero values are
    then quantile-binned among just the *other* nonzero values into
    classes `2..n_classes`, so the zero-heavy mass doesn't compress every
    real nonzero value into one or two top buckets."""
    if value <= 0:
        return 1
    nonzero = sorted(v for v in all_values if v > 0)
    if not nonzero:
        return n_classes  # value is the only nonzero one on record
    return _quantile_class(value, nonzero, n_classes=n_classes - 1) + 1


def classify_range(value: float, all_values: Sequence[float], *, n_classes: int = N_CLASSES) -> int:
    """Quantile-binned class (`1..n_classes`) of `value` among
    `all_values` -- for a metric that isn't zero-inflated (e.g. %identity,
    which is bounded and roughly continuous across a hit set)."""
    return _quantile_class(value, sorted(all_values), n_classes=n_classes)


def _quantile_class(value: float, sorted_values: Sequence[float], *, n_classes: int) -> int:
    n = len(sorted_values)
    if n <= 1:
        return n_classes
    # rank = fraction of sorted_values <= value (ties all land in the same
    # class). ceil, not floor/int -- e.g. with n_classes=5, a value whose
    # rank is anywhere in (0.2, 0.4] belongs in class 2, not class 1;
    # truncating would push it down a class right after every boundary.
    rank = sum(1 for v in sorted_values if v <= value) / n
    cls = math.ceil(rank * n_classes)
    return max(1, min(n_classes, cls))


def _family_segments(family: Optional[str]) -> List[str]:
    if not family or family == "(none reported)":
        return []
    return [s.strip() for s in family.split(".") if s.strip()]


def classify_family(family: Optional[str], seed_family: Optional[str], *, n_classes: int = N_CLASSES) -> int:
    """How deep `family`'s UniProt SIMILARITY hierarchy
    (superfamily -> family -> subfamily, e.g. "Belongs to the protein
    kinase superfamily. CMGC Ser/Thr protein kinase family. CDC2/CDKX
    subfamily") matches `seed_family`'s, from the top: exact match all the
    way to the seed's most specific level -> `n_classes` (e.g. same
    CDC2/CDKX subfamily as CDK20 itself); one level short (same family,
    different subfamily, e.g. MAP kinase vs. CDC2/CDKX) -> `n_classes - 1`;
    only the top-level superfamily matches -> the middle class; no overlap
    or missing data -> class 1."""
    seed_segs = _family_segments(seed_family)
    segs = _family_segments(family)
    if not segs or not seed_segs:
        return 1
    matched = 0
    for a, b in zip(segs, seed_segs):
        if a != b:
            break
        matched += 1
    if matched == 0:
        return 1
    if matched >= len(seed_segs):
        return n_classes
    if matched == len(seed_segs) - 1:
        return n_classes - 1
    if matched == 1:
        return (n_classes // 2) + 1
    return max(2, n_classes - (len(seed_segs) - matched))


@dataclass
class RankedHit:
    accession: str
    gene: str
    pct_identity: float
    n_templates: Optional[int]
    n_activities: Optional[int]
    family: str
    identity_class: int
    templates_class: int
    activity_class: int
    family_class: int
    max_score: int  # the highest score any row in this ranking could reach (n_classes^2 * count_classes^2) -- normalized_score's denominator

    @property
    def score(self) -> int:
        return self.identity_class * self.templates_class * self.activity_class * self.family_class

    @property
    def normalized_score(self) -> float:
        """`score` rescaled to 0.0-1.0 (`score / max_score`) -- lets scores
        be compared across rankings run with different `n_classes`/
        `count_classes`, where the raw integer score's own range differs."""
        return self.score / self.max_score


def _template_count(row: dict) -> Optional[int]:
    """The best RCSB-template-count signal available for a row: the exact,
    resolution-filtered count from an actual `--fetch`/`--fetch-all`
    (`pdb_structures`) if one has happened, else the cheaper, resolution-
    unfiltered total from `--pdb-count`/`--pdb-count-all` (`pdb_count`) if
    *that* has happened, else `None` ("not checked either way yet")."""
    if row["pdb_structures"] is not None:
        return len(row["pdb_structures"])
    return row.get("pdb_count")


def rank_hits(
    out_dir: Union[str, Path], *, n_classes: int = N_CLASSES, count_classes: int = N_COUNT_CLASSES,
) -> List[RankedHit]:
    """Rank every `role == "blast_hit"` row in `{out_dir}/hits.json` by
    the product of four ordinal scores (identity, RCSB template count,
    ChEMBL activity count, family relatedness to the seed), highest first.
    Identity and family use `n_classes` classes (1-5 by default); template
    and activity counts -- wider-ranging, less discrete signals -- use the
    finer `count_classes` (1-20 by default, see module docstring for why).
    RCSB template count prefers an exact, resolution-filtered `--fetch`/
    `--fetch-all` count when available, but falls back to the cheap,
    resolution-unfiltered `--pdb-count`/`--pdb-count-all` total so real
    template availability can inform the ranking *before* any actual
    structure download happens (see `_template_count`) -- the whole point
    of ranking ahead of fetching. Rows where neither RCSB signal nor
    `--chembl-activity` have been run yet (`pdb_structures`/`pdb_count`/
    `chembl_targets` still `None`) are scored as if that count were 0
    (class 1) -- "not yet checked" and "checked, found none" are
    deliberately not distinguished here, unlike in `hits.json` itself,
    since a ranking has to put every row somewhere."""
    out_dir = Path(out_dir)
    result = _load_hits_json(out_dir)
    seed = next(r for r in result["hits"] if r["role"] == "seed")
    hits = [r for r in result["hits"] if r["role"] == "blast_hit"]

    identities = [r["pct_identity"] for r in hits]
    template_counts = [_template_count(r) or 0 for r in hits]
    activity_counts = [sum(t["n_activities"] for t in (r["chembl_targets"] or [])) for r in hits]
    max_score = n_classes * count_classes * count_classes * n_classes

    ranked = []
    for r, n_templ, n_act in zip(hits, template_counts, activity_counts):
        ranked.append(RankedHit(
            accession=r["accession"], gene=r["gene"], pct_identity=r["pct_identity"],
            n_templates=_template_count(r),
            n_activities=None if r["chembl_targets"] is None else n_act,
            family=r["family"],
            identity_class=classify_range(r["pct_identity"], identities, n_classes=n_classes),
            templates_class=classify_zero_inflated(n_templ, template_counts, n_classes=count_classes),
            activity_class=classify_zero_inflated(n_act, activity_counts, n_classes=count_classes),
            family_class=classify_family(r["family"], seed["family"], n_classes=n_classes),
            max_score=max_score,
        ))
    ranked.sort(key=lambda h: h.score, reverse=True)
    return ranked
