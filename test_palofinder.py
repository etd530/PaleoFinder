#### LIBS ####
import paleofinder as pf
import filecmp

#### Tests for subset_gff ####
# Test that the right entries are kept from the unfiltered gff
def test_subset_gff():
	pf.os.system('mkdir -p paleofinder.extended_peptides_all_gff && cp pytest_inputs/paleofinder.OrNV_EU747721.1_ACH96265.1_135_protein_unknown.reconstructed_peptides.gff paleofinder.extended_peptides_all_gff')
	blastp_filtered_summary=pf.pd.read_csv("test_blastp_filtered_summary.tsv", sep='\t', header = 0)

	pf.subset_gff(blastp_filtered_summary=blastp_filtered_summary, outprefix="paleofinder")
	assert filecmp.cmp("pytest_inputs/test_gff_output.gff", "paleofinder.filtered_peptides_gff/paleofinder.OrNV_EU747721.1_ACH96265.1_135_protein_unknown.reconstructed_peptides.filtered.gff")
	pf.os.system("rm -rf extended_peptides_all_gff filtered_peptides_gff")

# ## Test for filtering of blastp output if original protein used to find the seeds is not in the blastp hits ##
# def filter_blastp_test_homolog_seqid():
# 	filter_blastp_output("blastp_test_output.out", 1, homologs_length_dict, os.path.dirname(__file__).strip('.') + '/nodes.dmp', os.path.dirname(__file__).strip('.') + '/names.dmp', os.path.dirname(__file__).strip('.') + '/merged.dmp', excluded_taxids_list = args['--excluded_taxids'])
