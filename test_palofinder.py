#!/usr/bin/env python3

#### LIBS ####
import paleofinder as pf
import filecmp
from Bio.Seq import Seq

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

#### Tests for extract_fragments_from_fasta ####
# Test that the correct fragments are extracted from the correct scaffolds
def test_extract_fragments_from_fasta():
	input_fragments_list = {'scaffold_1': [[[33, 43, 640, 719]]], 'scaffold_2': [[[38, 49, 126, 336]]], 'scaffold_3': [[[32, 50, 72, 161]]]}
	output_fragments_list = {'scaffold_1': [[[33, 43, Seq('GTAGCATGTTA'), 640, 719]]], 'scaffold_2': [[[38, 49, Seq('AGCATCTTATCG'), 126, 336]]], 'scaffold_3': [[[32, 50, Seq('CGACTAGCCGGCATAGGCG'), 72, 161]]]} 
	file = 'pytest_inputs/test_genome_input.fasta'

	pf.extract_fragments_from_fasta(file, input_fragments_list)
	assert input_fragments_list == output_fragments_list

#### Tests for extract_fragments_from_scaffold ####
# Test that the correct fragment is recovered when extracting from the forward strand
def test_extract_fragments_from_scaffold_forward_strand():
	coordinates = [15, 40]
	scaffold = Seq('TTACACAAAGGCATTTATGAGCGCGCTAGGGCTCACGAGCATCTTATCGAGCAGCGAGTTTAGCAGCATCACGAGCTACGTAGCGCGGTTTAGCAGCGAT')
	correct_fragment = Seq('TTATGAGCGCGCTAGGGCTCACGAGC')
	fragment = pf.extract_fragments_from_scaffold(scaffold, coordinates)
	assert fragment == correct_fragment

# Test that the correct fragment is recovered when extracting from the reverse strand
def test_extract_fragments_from_scaffold_reverse_strand():
	coordinates = [40, 15]
	scaffold = Seq('TTACACAAAGGCATTTATGAGCGCGCTAGGGCTCACGAGCATCTTATCGAGCAGCGAGTTTAGCAGCATCACGAGCTACGTAGCGCGGTTTAGCAGCGAT')
	correct_fragment = Seq('GCTCGTGAGCCCTAGCGCGCTCATAA')
	fragment = pf.extract_fragments_from_scaffold(scaffold, coordinates)
	assert fragment == correct_fragment

#### Tests for translate_dna ####
# Test that the correct translations are generated when using the first reading frame
def test_translate_dna_first_frame():
	input_fragments_list = {'scaffold_1': [[[33, 43, Seq('GTAGCATGTTA'), 640, 719]]], 'scaffold_2': [[[38, 49, Seq('AGCATCTTATCG'), 126, 336]]], 'scaffold_3': [[[32, 50, Seq('CGACTAGCCGGCATAGGCG'), 72, 161]]]} 
	correct_fragments_list = {'scaffold_1': [[[33, 43, Seq('GTAGCATGTTA'), Seq('VAC'), 640, 719]]], 'scaffold_2': [[[38, 49, Seq('AGCATCTTATCG'), Seq('SILS'), 126, 336]]], 'scaffold_3': [[[32, 50, Seq('CGACTAGCCGGCATAGGCG'), Seq('RLAGIG'), 72, 161]]]} 
	pf.translate_dna(input_fragments_list, frame=0)
	print(input_fragments_list)
	assert input_fragments_list == correct_fragments_list

#### Tests for align_peptides_simple ####
# NOT PRIORITARY TO TEST

#### Tests for align_bridges ####
# NOT PRIORITARY TO TEST

#### Tests for get_scaffold_from_fasta #####
def test_get_scaffold_from_fasta():
	correct_output = Seq('TACGTACGATAGCGCCTTAGCATGCATCGAGGGTAGCATGTTAGCGACGATCGATCAGCGGCGATCGATCGGACGATTCAGCATCGATTTTGCTAGCATG')
	output = pf.get_scaffold_from_fasta('pytest_inputs/test_genome_input.fasta', 'scaffold_1')
	assert output == correct_output

#### Tests for get_next_fragment_coordinates ####
# Test for correct coordinates when moving downstream and peptides being forward
def test_get_next_fragment_coordinates_forward_downstream():
	candidate_peptide = [10, 20, Seq('ATCGATCGATC'), Seq('MID'), 30, 32]
	direction = 'downstream'
	fragment_size = 30
	correct_start = 21
	correct_end = 50
	start, end = pf.get_next_fragment_coordinates(direction, candidate_peptide, fragment_size)
	assert start == correct_start and end == correct_end

# Test for correct coordinates when moving upstream and peptides being forward
def test_get_next_fragment_coordinates_forward_upstream():
	candidate_peptide = [100, 120, Seq('ATCGATCGATC'), Seq('MID'), 30, 32]
	direction = 'upstream'
	fragment_size = 30
	correct_start = 70
	correct_end = 99
	start, end = pf.get_next_fragment_coordinates(direction, candidate_peptide, fragment_size)
	assert start == correct_start and end == correct_end

####################################################################################
# Test for correct coordinates when moving downstream and peptides being reverse####
#### ADD TEST																	####
# Test for correct coordinates when moving upstream and peptides being reverse  ####
#### ADD TEST																	####
####################################################################################

#### Tests for check_fragment_contiguity ####
# def test_check_fragment_contiguity_upstream_lead_contiguous():
# 	input_candidate_peptide = [10, 20, Seq('ATCGATCGATC'), Seq('MID'), 30, 32]
# 	input_homolog_start = 


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
