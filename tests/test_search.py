"""Fast, offline unit tests for search.py's pure-logic helpers -- no
network access."""
from dd_idea.search import _gene_subdir


def test_gene_subdir_uses_gene_name():
    row = {"gene": "CDK2"}
    assert _gene_subdir("P24941", row) == "CDK2"


def test_gene_subdir_falls_back_to_accession_when_gene_unknown():
    row = {"gene": "-"}
    assert _gene_subdir("Q8IZL9", row) == "Q8IZL9"


def test_gene_subdir_falls_back_to_accession_when_gene_missing():
    row = {}
    assert _gene_subdir("Q8IZL9", row) == "Q8IZL9"
