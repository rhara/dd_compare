"""Fast, offline unit tests for rcsb/chain.py's chain-alignment/picking
logic -- ported alongside the vendored code they cover (from
dd_seqalign.sequence's own test_sequence.py -- see rcsb/chain.py's module
docstring), following this project family's standing pattern of porting a
vendored module's test coverage, not just its code."""
from dd_idea.rcsb import ChainSequence, align_to_canonical, pick_target_chain

CANONICAL = "MEDYTKIEKIGEGTYGVVYKGRHKTTGQVVAMKKIRLESEEEGVPSTAIREISLLKELRHPNIVSLQDVLMQDSRLYLIFEFLSMDLKKYLDSI"


def test_align_to_canonical_full_coverage_no_mismatch():
    chain = ChainSequence(chain_id="A", residues=list(enumerate(CANONICAL, start=1)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_covered == len(CANONICAL)
    assert aln.n_mismatch == 0
    assert aln.resseq_for_canonical(1) == 1


def test_align_to_canonical_fragment_has_free_end_gaps():
    fragment = CANONICAL[10:40]  # a co-crystal that only resolved residues 11-40
    chain = ChainSequence(chain_id="A", residues=list(enumerate(fragment, start=11)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_covered == 30
    assert aln.resseq_for_canonical(1) is None  # outside the fragment: missing, not a bad alignment
    assert aln.resseq_for_canonical(11) == 11
    assert aln.resseq_for_canonical(40) == 40


def test_align_to_canonical_flags_point_mutation_not_missing():
    mutant = CANONICAL[:20] + "A" + CANONICAL[21:]
    chain = ChainSequence(chain_id="A", residues=list(enumerate(mutant, start=1)))
    aln = align_to_canonical(chain, CANONICAL)
    assert aln.n_mismatch == 1
    mismatch = next(r for r in aln.residues if r.status == "mismatch")
    assert mismatch.canonical_pos == 21


def test_pick_target_chain_prefers_identity_over_raw_coverage():
    # A bound partner chain (e.g. a cyclin, a CAK) can align across the same
    # span as the true target chain while being mostly mismatches --
    # pick_target_chain must rank by matching residues, not raw coverage.
    true_chain = ChainSequence(chain_id="A", residues=list(enumerate(CANONICAL, start=1)))
    partner_seq = "".join("X" if i % 2 == 0 else c for i, c in enumerate(CANONICAL))
    partner_chain = ChainSequence(chain_id="B", residues=list(enumerate(partner_seq, start=1)))

    alignments = {
        "A": align_to_canonical(true_chain, CANONICAL),
        "B": align_to_canonical(partner_chain, CANONICAL),
    }
    assert pick_target_chain(alignments) == "A"
