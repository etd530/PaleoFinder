#!/usr/bin/env python3

#### LIBS ####
import filter_overlapping_features as fof

#### Tests for filter_overlapping_features ####
# Test that peptides cannot be discarded based on a previously discarded one
def test_filter_overlapping_features_discarding_by_discarded():
	gff = 'pytest_inputs/test_filter_overlapping_features_discarding_by_discarded.gff'
	blastp_summary = 'pytest_inputs/test_filter_overlapping_features_discarding_by_discarded.tsv'
	correct_peptides2keep = {'scaffold1___protein1___pseudopeptide_candidate_0', 'scaffold1___protein1___pseudopeptide_candidate_2'}
	peptides2keep = fof.filter_overlapping_features(gff, blastp_summary)
	print(peptides2keep)
	assert peptides2keep == correct_peptides2keep
