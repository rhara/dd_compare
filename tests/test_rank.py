"""Fast, offline unit tests for rank.py's classification/scoring logic --
no network access (rank_hits itself reads a pre-built hits.json fixture
via tmp_path, no live search.search()/fetch_templates() call)."""
import json

from dd_idea.rank import _template_count, classify_family, classify_range, classify_zero_inflated, rank_hits

CDK20_FAMILY = "Belongs to the protein kinase superfamily. CMGC Ser/Thr protein kinase family. CDC2/CDKX subfamily"
CDK2_FAMILY = CDK20_FAMILY  # same subfamily as CDK20
MAPK1_FAMILY = "Belongs to the protein kinase superfamily. CMGC Ser/Thr protein kinase family. MAP kinase subfamily"
PAK1_FAMILY = "Belongs to the protein kinase superfamily. STE Ser/Thr protein kinase family. STE20 subfamily"
NO_SUBFAMILY = "Belongs to the protein kinase superfamily. CMGC Ser/Thr protein kinase family"  # same family as CDK20, no subfamily on record


def test_classify_zero_inflated_zero_is_always_class_1():
    assert classify_zero_inflated(0, [0, 0, 5, 10, 100]) == 1


def test_classify_zero_inflated_nonzero_quantiles_among_nonzero_only():
    # Half the population is zero; nonzero values [5, 10, 100] shouldn't be
    # compressed by the zero mass -- the smallest nonzero still lands above
    # class 1, the largest lands at the top class.
    all_values = [0, 0, 0, 0, 0, 5, 10, 100]
    assert classify_zero_inflated(5, all_values, n_classes=5) >= 2
    assert classify_zero_inflated(100, all_values, n_classes=5) == 5


def test_classify_range_spans_full_class_range():
    values = [10, 20, 30, 40, 50, 60, 70, 80, 90, 100]
    assert classify_range(10, values, n_classes=5) == 1
    assert classify_range(100, values, n_classes=5) == 5
    # a value near the middle should land in a middle class, not get
    # truncated down to a lower one (regression check for the ceil-vs-int
    # rounding bug this function was first written with)
    assert classify_range(55, values, n_classes=5) in (3, 4)


def test_classify_range_single_value_gets_top_class():
    assert classify_range(42, [42], n_classes=5) == 5


def test_classify_family_exact_subfamily_match_is_top_class():
    assert classify_family(CDK2_FAMILY, CDK20_FAMILY, n_classes=5) == 5


def test_classify_family_same_family_different_subfamily_is_second_class():
    assert classify_family(MAPK1_FAMILY, CDK20_FAMILY, n_classes=5) == 4


def test_classify_family_same_superfamily_only_is_middle_class():
    assert classify_family(PAK1_FAMILY, CDK20_FAMILY, n_classes=5) == 3


def test_classify_family_missing_data_is_bottom_class():
    assert classify_family(None, CDK20_FAMILY) == 1
    assert classify_family("(none reported)", CDK20_FAMILY) == 1
    assert classify_family(CDK20_FAMILY, None) == 1


def test_classify_family_no_subfamily_on_record_still_credits_family_match():
    # NO_SUBFAMILY matches CDK20's first two segments (superfamily, family)
    # but has no third segment to compare -- one level short of an exact
    # subfamily match, same as MAPK1's case above.
    assert classify_family(NO_SUBFAMILY, CDK20_FAMILY, n_classes=5) == 4


def _hits_json(tmp_path, hits):
    out_dir = tmp_path
    (out_dir / "hits.json").write_text(json.dumps({
        "seed_accession": "SEED", "hits": [
            {"accession": "SEED", "gene": "SEED", "family": CDK20_FAMILY, "pct_identity": 100.0,
             "evalue": None, "pdb_structures": None, "chembl_targets": None, "role": "seed"},
            *hits,
        ],
    }))
    return out_dir


def _hit_row(acc, gene, family, pct_identity, n_templates, n_activities, *, pdb_count=None):
    return {
        "accession": acc, "gene": gene, "family": family, "pct_identity": pct_identity, "evalue": 1e-50,
        "pdb_structures": None if n_templates is None else [{"pdb_id": "X"}] * n_templates,
        "pdb_count": pdb_count,
        "chembl_targets": None if n_activities is None else [{"target_chembl_id": "T1", "n_activities": n_activities}],
        "role": "blast_hit",
    }


def test_rank_hits_orders_best_hit_first(tmp_path):
    out_dir = _hits_json(tmp_path, [
        _hit_row("WEAK", "WEAK1", PAK1_FAMILY, 25.0, 0, 0),
        _hit_row("STRONG", "CDK2", CDK2_FAMILY, 44.0, 275, 3015),
        _hit_row("MID", "MAPK1", MAPK1_FAMILY, 35.0, 5, 200),
    ])
    ranked = rank_hits(out_dir)
    assert [h.accession for h in ranked] == ["STRONG", "MID", "WEAK"]
    assert ranked[0].score > ranked[1].score > ranked[2].score


def test_normalized_score_is_score_over_max_possible(tmp_path):
    out_dir = _hits_json(tmp_path, [_hit_row("HIT1", "GENE1", CDK2_FAMILY, 44.0, 275, 3015)])
    ranked = rank_hits(out_dir, n_classes=5, count_classes=20)
    assert ranked[0].max_score == 5 * 20 * 20 * 5
    assert ranked[0].normalized_score == ranked[0].score / (5 * 20 * 20 * 5)


def test_normalized_score_top_hit_is_1_0(tmp_path):
    # A hit that's the sole nonzero value for every zero-inflated signal and
    # the seed's exact subfamily should max out every class -> normalized 1.0.
    out_dir = _hits_json(tmp_path, [_hit_row("HIT1", "GENE1", CDK2_FAMILY, 44.0, 1, 1)])
    ranked = rank_hits(out_dir)
    assert ranked[0].normalized_score == 1.0


def test_normalized_score_comparable_across_different_count_classes(tmp_path):
    # The whole point of normalizing: two rankings run with different
    # count_classes settings still land the same best-possible hit at 1.0.
    out_dir = _hits_json(tmp_path, [_hit_row("HIT1", "GENE1", CDK2_FAMILY, 44.0, 1, 1)])
    ranked_5 = rank_hits(out_dir, count_classes=5)
    ranked_20 = rank_hits(out_dir, count_classes=20)
    assert ranked_5[0].score != ranked_20[0].score  # raw scores differ...
    assert ranked_5[0].normalized_score == ranked_20[0].normalized_score == 1.0  # ...but normalized agree


def test_rank_hits_treats_unfetched_rows_as_zero_not_missing(tmp_path):
    out_dir = _hits_json(tmp_path, [
        _hit_row("NOTFETCHED", "X1", CDK2_FAMILY, 40.0, None, None),
    ])
    ranked = rank_hits(out_dir)
    assert ranked[0].n_templates is None
    assert ranked[0].n_activities is None
    assert ranked[0].templates_class == 1
    assert ranked[0].activity_class == 1


def test_rank_hits_excludes_seed_row(tmp_path):
    out_dir = _hits_json(tmp_path, [_hit_row("HIT1", "GENE1", CDK2_FAMILY, 40.0, 1, 1)])
    ranked = rank_hits(out_dir)
    assert "SEED" not in [h.accession for h in ranked]


def test_rank_hits_default_count_classes_separate_ties_5_classes_would_collapse(tmp_path):
    # CDK2's real numbers (522 templates/3015 activities) and CDK9's (28/2051)
    # both landed in template/activity class 5 under the old 5-class scheme,
    # tying their scores despite being nowhere close -- the whole motivation
    # for widening templates_class/activity_class's default granularity.
    out_dir = _hits_json(tmp_path, [
        _hit_row("CDK2", "CDK2", CDK2_FAMILY, 43.8, 522, 3015),
        _hit_row("CDK9", "CDK9", CDK2_FAMILY, 40.0, 28, 2051),
        *[_hit_row(f"FILLER{i}", f"F{i}", PAK1_FAMILY, 25.0, i, i * 10) for i in range(1, 10)],
    ])
    ranked = rank_hits(out_dir)
    cdk2 = next(h for h in ranked if h.accession == "CDK2")
    cdk9 = next(h for h in ranked if h.accession == "CDK9")
    assert cdk2.templates_class > cdk9.templates_class
    assert cdk2.activity_class > cdk9.activity_class
    assert cdk2.score > cdk9.score


def test_rank_hits_count_classes_is_overridable(tmp_path):
    out_dir = _hits_json(tmp_path, [
        _hit_row("A", "A1", CDK2_FAMILY, 40.0, 10, 10),
        _hit_row("B", "B1", CDK2_FAMILY, 40.0, 20, 20),
    ])
    ranked = rank_hits(out_dir, count_classes=5)
    assert all(h.templates_class <= 5 and h.activity_class <= 5 for h in ranked)


def test_template_count_prefers_fetched_over_pdb_count():
    # pdb_structures (an actual --fetch, resolution-filtered) wins over a
    # cheaper --pdb-count total even when both are on record.
    row = _hit_row("X", "X1", CDK2_FAMILY, 40.0, n_templates=2, n_activities=None, pdb_count=500)
    assert _template_count(row) == 2


def test_template_count_falls_back_to_pdb_count_when_not_fetched():
    row = _hit_row("X", "X1", CDK2_FAMILY, 40.0, n_templates=None, n_activities=None, pdb_count=500)
    assert _template_count(row) == 500


def test_template_count_none_when_neither_checked():
    row = _hit_row("X", "X1", CDK2_FAMILY, 40.0, n_templates=None, n_activities=None, pdb_count=None)
    assert _template_count(row) is None


def test_rank_hits_uses_pdb_count_before_any_fetch_has_happened(tmp_path):
    # This is the whole point of --pdb-count: a hit set can be meaningfully
    # ranked on real RCSB template availability *before* --fetch/--fetch-all
    # ever runs, unlike when pdb_structures alone (still None for everyone)
    # would leave every row tied at templates_class 1.
    out_dir = _hits_json(tmp_path, [
        _hit_row("MANY", "CDK2", CDK2_FAMILY, 44.0, n_templates=None, n_activities=None, pdb_count=500),
        _hit_row("NONE", "MAK", CDK2_FAMILY, 35.0, n_templates=None, n_activities=None, pdb_count=0),
    ])
    ranked = rank_hits(out_dir)
    many = next(h for h in ranked if h.accession == "MANY")
    none = next(h for h in ranked if h.accession == "NONE")
    assert many.n_templates == 500
    assert none.n_templates == 0
    assert many.templates_class > none.templates_class
    assert many.score > none.score
