#!/usr/bin/env python3

"""
Filter the output of PaleoFinder based on the overlap among the reconstructed peptides.

Usage:
	filter_overlapping_features.py --gff=GFF --fasta=FASTA --blastp=BLASTP --blastp_summary=TSV
	filter_overlapping_features.py -h | --help

Arguments:
	--gff                                             Filtered GFF file ouput from PaleoFinder to filter out.
	--fasta                                           FASTA file to filter.
	--blastp                                          BLASTP tabular output file to filter.
	--blastp_summary                                  BLASTP summary output file to filter.

Options:
	-h --help:                                        Show this message and exit.
"""

#### LIBS ####
from filter_repetitive_proteins import subset_blastp_output, subset_blastp_summary, subset_fasta, subset_gff
from docopt import docopt
import re
import portion as P

#### FUNS ####
def filter_overlapping_features(gff, blastp_summary):
	"""
	Given a GFF file, filter those features with overlapping coordinates, keeping the one with the highest Alien Index form the TSV summary file.

	Arguments:
		gff: GFF file to filter.
		blastp_summary: TSV file with Alien Index values.
	
	Returns:
		peptides2keep: a list with peptide IDs that need to be kept.
	"""
	def get_peptide_name(line):
		gff_entry = line.strip().split('\t')
		scaffold = gff_entry[0]
		homolog = gff_entry[-1].split(';')[0].replace('ID=', '')
		peptide_number = re.sub(r".+-like\.pseudogene_([0-9]+)\.[0-9]+", r"pseudopeptide_candidate_\1", homolog)
		homolog = re.sub(r"-like\.pseudogene_([0-9]+)\.[0-9]+", '', homolog)
		peptide_name = '___'.join([scaffold, homolog, peptide_number])
		return peptide_name

	# Get a dictionary with Alien Index values of all peptides
	alien_index_dict = {}
	with open(blastp_summary, 'r') as fh:
		for line in fh:
			if line.startswith('query'):
				continue
			else:
				peptide = line.strip().split('\t')[0]
				alien_index = float(line.strip().split('\t')[6])
				alien_index_dict[peptide] = alien_index
	
	# Compare the peptides in the GFF file and make a list with non-overlapping peptides, choosing the one with highest Alien Index.
	peptides_set = set() # to store all peptides
	peptides2remove = set() # to store those that overlap at some point
	
	with open(gff, 'r') as fh:
		for line in fh:
			if line.startswith('##'):
				continue
			else:
				discard = False # a priori each line should be kept
				peptide_name1 = get_peptide_name(line)
				if peptide_name1 in peptides2remove: # if the peptide has been discarded, skip it
					continue
				peptides_set.add(peptide_name1)
				alien_index1 = alien_index_dict[peptide_name1]
				range1 = [int(line.strip().split('\t')[3]), int(line.strip().split('\t')[4])]
				len1 = len(range(range1[0], range1[1]))
				scaffold1 = line.strip().split('\t')[0]
				with open(gff, 'r') as fh2:
					for line2 in fh2:
						if line2.startswith('##'):
							continue
						else:
							peptide_name2 = get_peptide_name(line2)
							if peptide_name1 == peptide_name2 or peptide_name2 in peptides2remove: # skip comparisons of fragments of the same peptide or involving peptides already discarded
								continue
							else:
								scaffold2 = line2.strip().split('\t')[0]
								if scaffold1 == scaffold2: # if both peptides are on the same scaffold check coordinates overlap
									range2 = [int(line2.strip().split('\t')[3]), int(line2.strip().split('\t')[4])]
									overlap = P.closed(range1[0], range1[1]).intersection(P.closed(range2[0], range2[1]))
									if overlap: # if there is overlap check the alien index
										alien_index2 = alien_index_dict[peptide_name2]
										if alien_index1 < alien_index2: # if alien index is lower, discard this peptide
											peptides2remove.add(peptide_name1)
											print("Discarding peptide %s because Alien Index is lower than that of %s." % (peptide_name1, peptide_name2))
											discard = True
										elif alien_index1 == alien_index2: # if both AIs are equal, compare the lengths
											len2 = len(range(range2[0], range2[1]))
											if len1 < len2: # if the peptide is shorter, discard it
												peptides2remove.add(peptide_name1)
												print("Discarding peptide %s because it is shorter than %s" % (peptide_name1, peptide_name2))
												discard = True
											elif len1 == len2 and peptide_name2 not in peptides2remove:
												# if both have the same length and the 2nd one has not been discarded, discard also the first one
												peptides2remove.add(peptide_name1)
												print("Tie. Discarding peptide %s because the other one (%s) has not been discarded yet." % (peptide_name1, peptide_name2))
												discard = True
							if discard: # if the peptide is discarded, break the loop
								break				
	# After determining the peptides to be discarded, the ones to be kept is all except those in the discard set
	peptides2keep = peptides_set - peptides2remove
	return peptides2keep

if __name__ == '__main__':
	#### ARGUMENTS ####
	args = docopt(__doc__)
	gff = args['--gff']
	fasta = args['--fasta']
	blastp = args['--blastp']
	blastp_summary = args['--blastp_summary']

	#### MAIN ####
	# Get list of non-overlapping peptides to be kept
	peptides2keep = filter_overlapping_features(gff, blastp_summary)
	print(peptides2keep)

	subset_gff(gff, peptides2keep, '.filtered4overlapping.gff')
	subset_fasta(fasta, peptides2keep, '.filtered4overlapping.fasta')
	subset_blastp_output(blastp, peptides2keep, '.filtered4overlapping.out')
	subset_blastp_summary(blastp_summary, peptides2keep, '.filtered4overlapping.tsv')
