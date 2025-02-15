#!/usr/bin/env python3

#### LIBS ####
import filter_repetitive_proteins as frp

#### Tests for gff_ids_list ####
# Test that all the protein IDs are extracted from the GFF file
def test_list_peptides():
	correct_output = ['scaffold1___protein1___pseudopeptide_candidate_1', 'scaffold3___protein3___pseudopeptide_candidate_0', 'scaffold1___protein2___pseudopeptide_candidate_1', 'scaffold2___protein3___pseudopeptide_candidate_2']
	output = frp.list_peptides('test_fasta_files.fasta')
	print(output)
	assert output == correct_output

#### Tests for filter_peptides ####
# Test that the correct list of peptides to be kept is returned
def test_filter_peptides():
	correct_output = ['scaffold1___protein2___pseudopeptide_candidate_1', 'scaffold2___protein3___pseudopeptide_candidate_2', 'scaffold3___protein3___pseudopeptide_candidate_0']
	gff_ids_list = ['scaffold1___protein1___pseudopeptide_candidate_1', 'scaffold1___protein2___pseudopeptide_candidate_1', 'scaffold2___protein3___pseudopeptide_candidate_2', 'scaffold3___protein3___pseudopeptide_candidate_0']
	output = frp.filter_peptides(gff_ids_list, 'test_interproscan_file.gff3')

	# Confirm that the protien with no annotations is in the output
	assert 'scaffold3___protein3___pseudopeptide_candidate_0' in output

	# Confirm all the rest is also OK
	assert output == correct_output
