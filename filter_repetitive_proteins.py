#!/usr/bin/env python3

"""
Filter the output of PaleoFinder based on the presence of ankyrin repeats.

Usage:
	filter_repetitive_proteins.py --interproscan=GFF --gff=GFF --fasta=FASTA --blastp=BLASTP --blastp_summary=TSV
	filter_repetitive_proteins.py -h | --help

Arguments:
	--interproscan GFF                                GFF file with Interproscan annotations.
	--gff                                             Filtered GFF file ouput from PaleoFinder to filter out.
	--fasta                                           FASTA file to filter.
	--blastp                                          BLASTP tabular output file to filter.
	--blastp_summary                                  BLASTP summary output file to filter.

Options:
	-h --help:                                        Show this message and exit.
"""

#### LIBS ####
from docopt import docopt
from BCBio import GFF
import pprint
import portion as P
import re
import glob
from Bio.SeqIO import FastaIO
from Bio import SeqIO
import os

#### FUNCTIONS ####
def list_peptides(fasta):
	"""
	Given a FASTA file, extract the list of protein IDs in it.

	Arguments:
		fasta: FASTA file.
	
	Returns:
		gff_ids_list: List of protein IDs
	"""
	gff_ids_list = []
	# Open the FASTA file
	with open(fasta, "r") as fh:
		# Iterate over each sequence in the FASTA file
		for record in SeqIO.parse(fh, "fasta"):
			# Access the sequence ID
			gff_ids_list.append(record.id)
	return gff_ids_list

def remove_overlap(ankyrin_features):
	"""
	Remove overlap among the annotated ankyrin features.

	Arguments:
		ankyrin_features: List of lists with the start and end of the ankyrin features.

	Returns:
		final_range: List of lists with the non-overlapping ankyrin features.
	"""
	intervals = [P.closed(a, b) for a, b in ankyrin_features]
	# Merge these intervals
	merge = P.Interval(*intervals)
	# Output as a list of lists
	final_range = [[i.lower, i.upper] for i in merge]
	return final_range

def filter_peptides(gff_ids_list, interproscan_gff):
	"""
	Filter peptides based on the presence of ankyrin repeats consistint in more than 50% of the protein.

	Arguments:
		gff_ids_list: List of protein IDs to filter
		interproscan_gff: GFF file with Interproscan annotations.
	
	Returns:
		peptides2keep: List of protein IDs to keep
	"""
	peptides2keep = [] # list to store the peptides that will be kept
	for id in gff_ids_list:
		print(id)
		ankyrin_features = [] # to store the ranges of ankryin features
		limit_info = dict(gff_id=[id])
		in_handle = open(interproscan_gff)
		for rec in GFF.parse(in_handle, limit_info=limit_info):
			for feature in rec.features:
				if hasattr(feature, 'qualifiers'):
					if feature.type == 'polypeptide': # polypeptide just gives the whole peptide, extract the length
						peptide_start = feature.location.start
						peptide_end = feature.location.end
						peptide_len = len(range(peptide_start, peptide_end))
					elif feature.type == 'protein_match' and 'signature_desc' in feature.qualifiers:
						signature_desc = '_'.join(feature.qualifiers['signature_desc'])
						if 'Ankyrin' in signature_desc: # if Ankyrin is annotated, get the coordinates
							feature_start = feature.location.start
							feature_end = feature.location.end
							ankyrin_features.append([feature_start, feature_end])
		# Once all the features have been extracted, remove overlap among them and get the total length
		if len(ankyrin_features) > 1:
			print("Has ankyrin repeats")
			ankyrin_features = remove_overlap(ankyrin_features)
			print(ankyrin_features)
			ankyrin_len = sum([len(range(i[0], i[1])) for i in ankyrin_features])
			ankyrin_percentage = (ankyrin_len / peptide_len) * 100
			# If ankyrin percentage is lower than 50% we keep the protein
			if ankyrin_percentage < 50:
				print("Has less than 50% ankyrin repeats")
				peptides2keep.append(id)
			else:
				print("Has more than 50% ankyrin repeats")
		else:
			print("Does not have ankyrin repeats")
			peptides2keep.append(id)
		in_handle.close()
	return peptides2keep

def subset_gff(gff, peptides2keep, outsuffix):
	"""
	Write a new GFF file containing only those peptides that are in the input list.

	Arguments:
		gff: GFF file to filter.
		peptides2keep: List object containing the protein IDs to keep.
		outsuffix: Suffix to add to the output file.

	Returns:
		None. Writes the GFF files to disk.
	"""
	outprefix=os.path.basename(gff).strip('.gff3')
	with open(gff, 'r') as fh:
		for line in fh:
			# print(line)
			if line.startswith('##'):
				file_string = line
			else:
				gff_entry = line.strip().split('\t')
				scaffold = gff_entry[0]
				homolog = gff_entry[-1].split(';')[0].replace('ID=', '')
				peptide_number = re.sub(r".+-like\.pseudogene_([0-9]+)\.[0-9]+", r"pseudopeptide_candidate_\1", homolog)
				homolog = re.sub(r"-like\.pseudogene_([0-9]+)\.[0-9]+", '', homolog)
				peptide_name = '___'.join([scaffold, homolog, peptide_number])
				print(peptide_name)
				if peptide_name in peptides2keep:
					print('Peptide needs to be kept')
					file_string = file_string + line
				else:
					print('Peptide does not need to be kept')
	if file_string != '##gff-version 3\n':
		with open(outprefix+outsuffix, 'w') as wh:
			wh.write(file_string)

def subset_fasta(fasta, peptides2keep, outsuffix):
	"""
	Write a new set of FASTA files containing only those peptides that are in the input list.

	Arguments:
		fasta: Folder containing the FASTA files to filter.
		peptides2keep: List object containing the protein IDs to keep.
		outsuffix: Suffix to add to the output file.

	Returns:
		None. Writes the FASTA files to disk.
	"""
	file_string = ''
	with open(fasta, 'r') as fh:
		for sequence in FastaIO.FastaIterator(fh):
			if sequence.id in peptides2keep:
				file_string = file_string + ">%s\n%s\n" % (sequence.id, sequence.seq)
	if len(file_string) > 0:
		outfile = os.path.basename(fasta).replace('.fasta', outsuffix)
		with open(outfile, 'w') as wh:
			wh.write(file_string)
	else:
		print('Filtered FASTA file is empty!')

def subset_blastp_output(blastp, peptides2keep, outsuffix):
	"""
	Filter BLAST tabular output and summary output to keep only the peptides that are in the input list.

	Arguments:
		blastp: BLASTP tabular output file.
		peptides2keep: List object containing the protein IDs to keep.
		outsuffix: Suffix to add to the output file.
	
	Returns:
		None. Writes the filtered BLASTP files to disk.
	"""
	file_string = ''
	with open(blastp, 'r') as fh:
		for line in fh:
			if line.startswith('qseqid'):
				file_string = line
			else:
				peptide = line.strip().split('\t')[0]
				if peptide in peptides2keep:
					file_string = file_string + line
	with open(os.path.basename(blastp).replace('.out', outsuffix), 'w') as wh:
		wh.write(file_string)

def subset_blastp_summary(blastp_summary, peptides2keep, outsuffix):
	"""
	Filter TSV BLASTP summary files to keep only the peptides that are in the input list.

	Arguments:
		blastp_summary: BLASTP summary output file.
		peptides2keep: List object containing the protein IDs to keep.
		outsuffix: Suffix to add to the output file.

	Returns:
		None. Writes the filtered BLASTP summary file to disk.
	"""
	file_string = ''
	with open(blastp_summary, 'r') as fh:
		for line in fh:
			if line.startswith('query'):
				file_string = line
			else:
				peptide = line.strip().split('\t')[0]
				if peptide in peptides2keep:
					file_string = file_string + line
	with open(os.path.basename(blastp_summary).replace('.tsv', outsuffix), 'w') as wh:
		wh.write(file_string)

if __name__ == '__main__':
	#### ARGUMENTS ####
	args = docopt(__doc__)
	interproscan_gff = args['--interproscan']
	gff = args['--gff']
	fasta = args['--fasta']
	blastp = args['--blastp']
	blastp_summary = args['--blastp_summary']

	#### MAIN ####
	# Extract the list of peptides	
	gff_ids_list = list_peptides(fasta)

	# For each peptide extract the protein_match features and check if it is an Ankyrin
	# then check if it has too much Ankyirin and needs to be discarded
	peptides2keep = filter_peptides(gff_ids_list, interproscan_gff)
	print("LIST OF PEPTIDES TO KEEP:")
	print(peptides2keep)

	# Subset the GFF, FASTA, BLASTP and TSV summary files
	subset_gff(gff, peptides2keep, '.filtered4repeats.gff')
	subset_fasta(fasta, peptides2keep, '.filtered4repeats.fasta')
	subset_blastp_output(blastp, peptides2keep, '.filtered4repeats.out')
	subset_blastp_summary(blastp_summary, peptides2keep, '.filtered4repeats.tsv')
