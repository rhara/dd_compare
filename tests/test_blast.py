"""Fast, offline unit tests for blast.py's query classification and BLAST-
XML hit parsing -- no real network access (the XML fixture below stands in
for a real NCBI QBLAST response)."""
import pytest

from dd_idea import blast, search

_XML_FIXTURE = """<?xml version="1.0"?>
<!DOCTYPE BlastOutput PUBLIC "-//NCBI//NCBI BlastOutput/EN" "http://www.ncbi.nlm.nih.gov/dtd/NCBI_BlastOutput.dtd">
<BlastOutput>
  <BlastOutput_program>blastp</BlastOutput_program>
  <BlastOutput_version>BLASTP 2.13.0+</BlastOutput_version>
  <BlastOutput_db>swissprot</BlastOutput_db>
  <BlastOutput_query-ID>Query_1</BlastOutput_query-ID>
  <BlastOutput_query-def>unnamed protein product</BlastOutput_query-def>
  <BlastOutput_query-len>346</BlastOutput_query-len>
  <BlastOutput_param>
    <Parameters>
      <Parameters_matrix>BLOSUM62</Parameters_matrix>
      <Parameters_expect>1e-10</Parameters_expect>
      <Parameters_gap-open>11</Parameters_gap-open>
      <Parameters_gap-extend>1</Parameters_gap-extend>
      <Parameters_filter>F</Parameters_filter>
    </Parameters>
  </BlastOutput_param>
  <BlastOutput_iterations>
    <Iteration>
      <Iteration_iter-num>1</Iteration_iter-num>
      <Iteration_query-ID>Query_1</Iteration_query-ID>
      <Iteration_query-len>346</Iteration_query-len>
      <Iteration_hits>
        <Hit>
          <Hit_num>1</Hit_num>
          <Hit_id>sp|P24941.1|CDK2_HUMAN</Hit_id>
          <Hit_def>RecName: Full=Cyclin-dependent kinase 2</Hit_def>
          <Hit_accession>P24941</Hit_accession>
          <Hit_len>298</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>250</Hsp_bit-score>
              <Hsp_score>640</Hsp_score>
              <Hsp_evalue>1.1e-79</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>44</Hsp_identity>
              <Hsp_positive>60</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>ACDEFGHIKL</Hsp_qseq>
              <Hsp_hseq>ACDEFGHIKL</Hsp_hseq>
              <Hsp_midline>ACDEFGHIKL</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>2</Hit_num>
          <Hit_id>sp|P20794.2|MAK_HUMAN</Hit_id>
          <Hit_def>RecName: Full=Serine/threonine-protein kinase MAK</Hit_def>
          <Hit_accession>P20794</Hit_accession>
          <Hit_len>623</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>180</Hsp_bit-score>
              <Hsp_score>450</Hsp_score>
              <Hsp_evalue>4.2e-58</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>35</Hsp_identity>
              <Hsp_positive>50</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>ACDEFGHIKL</Hsp_qseq>
              <Hsp_hseq>ACDEFGHIKL</Hsp_hseq>
              <Hsp_midline>ACDEFGHIKL</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
        <Hit>
          <Hit_num>3</Hit_num>
          <Hit_id>sp|Q8IZL9.1|CDK20_HUMAN</Hit_id>
          <Hit_def>RecName: Full=Cyclin-dependent kinase 20 (self hit)</Hit_def>
          <Hit_accession>Q8IZL9</Hit_accession>
          <Hit_len>346</Hit_len>
          <Hit_hsps>
            <Hsp>
              <Hsp_num>1</Hsp_num>
              <Hsp_bit-score>700</Hsp_bit-score>
              <Hsp_score>1800</Hsp_score>
              <Hsp_evalue>0.0</Hsp_evalue>
              <Hsp_query-from>1</Hsp_query-from>
              <Hsp_query-to>100</Hsp_query-to>
              <Hsp_hit-from>1</Hsp_hit-from>
              <Hsp_hit-to>100</Hsp_hit-to>
              <Hsp_query-frame>0</Hsp_query-frame>
              <Hsp_hit-frame>0</Hsp_hit-frame>
              <Hsp_identity>100</Hsp_identity>
              <Hsp_positive>100</Hsp_positive>
              <Hsp_gaps>0</Hsp_gaps>
              <Hsp_align-len>100</Hsp_align-len>
              <Hsp_qseq>ACDEFGHIKL</Hsp_qseq>
              <Hsp_hseq>ACDEFGHIKL</Hsp_hseq>
              <Hsp_midline>ACDEFGHIKL</Hsp_midline>
            </Hsp>
          </Hit_hsps>
        </Hit>
      </Iteration_hits>
      <Iteration_stat>
        <Statistics>
          <Statistics_db-num>1</Statistics_db-num>
          <Statistics_db-len>1</Statistics_db-len>
          <Statistics_hsp-len>0</Statistics_hsp-len>
          <Statistics_eff-space>0</Statistics_eff-space>
          <Statistics_kappa>0.041</Statistics_kappa>
          <Statistics_lambda>0.267</Statistics_lambda>
          <Statistics_entropy>0.14</Statistics_entropy>
        </Statistics>
      </Iteration_stat>
    </Iteration>
  </BlastOutput_iterations>
</BlastOutput>
"""


def test_parse_blast_hits_excludes_seed_and_ranks_by_identity():
    hits = blast.parse_blast_hits(_XML_FIXTURE, exclude_accession="Q8IZL9", max_hits=10)
    accessions = [h.accession for h in hits]
    assert "Q8IZL9" not in accessions
    assert accessions == ["P24941", "P20794"]  # 44% before 35%
    assert hits[0].pct_identity == pytest.approx(44.0)
    assert hits[1].pct_identity == pytest.approx(35.0)


def test_parse_blast_hits_respects_max_hits():
    hits = blast.parse_blast_hits(_XML_FIXTURE, exclude_accession="Q8IZL9", max_hits=1)
    assert len(hits) == 1
    assert hits[0].accession == "P24941"


def test_parse_blast_hits_keeps_description_and_evalue():
    hits = blast.parse_blast_hits(_XML_FIXTURE, exclude_accession="Q8IZL9", max_hits=10)
    mak = next(h for h in hits if h.accession == "P20794")
    assert "MAK" in mak.description
    assert mak.evalue == pytest.approx(4.2e-58)


@pytest.mark.parametrize("query,expected", [
    ("Q8IZL9", "uniprot"),
    ("P24941", "uniprot"),
    ("q8izl9", "uniprot"),
    ("CHEMBL301", "chembl"),
    ("chembl301", "chembl"),
    ("MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQ", "sequence"),
])
def test_classify_query_recognizes_each_form(query, expected):
    assert search.classify_query(query) == expected


@pytest.mark.parametrize("query", ["", "short", "12345", "NOT_A_VALID_QUERY!!"])
def test_classify_query_rejects_garbage(query):
    with pytest.raises(ValueError):
        search.classify_query(query)
