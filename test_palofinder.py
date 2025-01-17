#!/usr/bin/env python3

#### LIBS ####
import paleofinder as pf
import filecmp

#### Tests for make_blast_db ####
# NOT PRIORITARY TO TEST

#### Tests for tblastn ####
# NOT PRIORITARY TO TEST

#### Tests for find_primary_seed_coordinates ####
# Test that the right primary seed coordinates are found
def test_find_primary_seed_coordinates():
	correct_seed_coordinates = {'scaffold_844': [[[33761, 33522, 640, 719]]], 'scaffold_872': [[[3852, 4493, 126, 336]]], 'scaffold_47': [[[318879, 318601, 72, 161]]]}
	
	input_tblastn_results = pf.pd.read_csv("pytest_inputs/test_tblastn_input.out", sep='\t', header = None)
	input_tblastn_results.rename(columns={0: 'qseqid', 1: 'sseqid', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'}, inplace = True)

	seed_coordinates = pf.find_primary_seed_coordinates(input_tblastn_results, query_seq_id='OrNV_EU747721.1_ACH96131.1_1_[protein=dnapol_B]')
	print(seed_coordinates)
	assert seed_coordinates == correct_seed_coordinates	

#### Tests for subset_gff ####
# Test that the right entries are kept from the unfiltered gff
def test_subset_gff():
	pf.os.system('mkdir -p paleofinder.extended_peptides_all_gff && cp pytest_inputs/test_gff_input.gff paleofinder.extended_peptides_all_gff/reconstructed_peptides.gff')
	blastp_filtered_summary=pf.pd.read_csv("pytest_inputs/test_blastp_filtered_summary.tsv", sep='\t', header = 0)

	pf.subset_gff(blastp_filtered_summary=blastp_filtered_summary, outprefix="paleofinder")
	assert filecmp.cmp("pytest_outputs/test_gff_output.gff", "paleofinder.filtered_peptides_gff/reconstructed_peptides_filtered.gff")
	pf.os.system("rm -rf paleofinder.extended_peptides_all_gff paleofinder.filtered_peptides_gff")

# ## Test for filtering of blastp output if original protein used to find the seeds is not in the blastp hits ##
# def filter_blastp_test_homolog_seqid():
# 	filter_blastp_output("blastp_test_output.out", 1, homologs_length_dict, os.path.dirname(__file__).strip('.') + '/nodes.dmp', os.path.dirname(__file__).strip('.') + '/names.dmp', os.path.dirname(__file__).strip('.') + '/merged.dmp', excluded_taxids_list = args['--excluded_taxids'])
