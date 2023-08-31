#!/usr/bin/env python3


"""
usage: pseudogene_finder.py --proteins=FASTA --genome=FASTA [--tblastn_wordsize=INT --tblastn_matrix=STR --tblastn_max_evalue=FLT --tblastn_seg_filter=STR --tblastn_threads=INT --outprefix=STR] [-v|--verbose] [-h|--help]

    Options:
        -p, --proteins FASTA                      Input query proteins for that serve as reference to find pseudogenes, in FASTA format. Can contain multiple sequences
        -g, --genome FASTA                        Input target genome in which pseudogenes are to be found, in FASTA format.
        --tblastn_wordsize INT                    Word size to use for the initial tBLASTN search [default: 3].
        --tblastn_matrix STR                      Alignment scoring matrix to use for the initial tBLASTn search [default: BLOSUM45].
        --tblastn_max_evalue FLT                  Maximum e-value to keep a hit in the initial tBLASTn search [default: 50].
        --tblastn_seg_filter STR                  Parameters for the SEG masking of low complexity regions in the query proteins [default: "10 1.0 1.5"].
        --tblastn_threads INT                     Number of threads to use to run the initial tBLASTN search [default: 10].
        --o, --outprefix STR                      Prefix to use for output files [default: pseudogene_finder].
        -v, --verbose                             Print the progressions of the program to the terminal (Standard Error).
        -h, --help                                Show this help message and exit.
"""


#### LIBS ####
from docopt import docopt
from Bio.SeqIO import FastaIO
from Bio import Align, AlignIO
import os
import pandas as pd


#### FUNS ####
def make_blast_db(genome):
	if not os.path.exists(genome+'.nhr'):
		os.system("makeblastdb -dbtype nucl -in " + genome)
	else:
		print("WARNING: BLAST database already found for " + genome + ". Skipping database building.") # Could implement a force flag if they want to rerun, but let's see if we end up using submodules...


def tblastn(query, target, wordsize, matrix, max_evalue, seg_filter, threads, outprefix):
	blast_file = outprefix + ".wordsize" + wordsize + "." + matrix + ".evalue" + max_evalue + ".seg" + seg_filter.replace('"', '') + ".out"
	blast_file = blast_file.replace(' ', '_')

	if not os.path.exists(blast_file):
		tblastn_command="tblastn -query " + query + " -db " + target + " -word_size " + wordsize + " -matrix " + matrix + " -evalue " + max_evalue + " -outfmt 6 -num_threads " + threads + " -seg " + seg_filter + " -out " + blast_file
		os.system(tblastn_command)
	else:
		print("WARNING: tBLASTn results found at " + blast_file + ". Skipping tBLASTN step.")
	
	tblastn_results = pd.read_csv(blast_file, sep='\t', header = None)
	tblastn_results.rename(columns={0: 'qseqid', 1: 'sseqid', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'}, inplace = True)
	
	return(tblastn_results)

def find_primary_seed_coordinates(tblastn_df, query_seq_id):
	tblastn_df_subset = tblastn_df[tblastn_df['qseqid'] == query_seq_id]
	primary_seeds = {}
	for index, row in tblastn_df_subset.iterrows():
		newkey = row['sseqid']
		if newkey in primary_seeds.keys():
			primary_seeds[newkey].append([[row['sstart'], row['send']]])
		else:
			primary_seeds[newkey] = [[[row['sstart'], row['send']]]]
		
		# print(newkey)
		# primary_seeds

	# = list(tblastn_df_subset.apply(lambda row: [row['sseqid'], row['sstart'], row['send']], axis=1))
	return(primary_seeds)

# def load_full_fasta(file):
# 	fasta_dict = {}
# 	with open(file) as fh:
# 		for sequence in FastaIO.FastaIterator(fh):
# 			fasta_dict[sequence.id] = sequence.seq
#  NOT FINISHED; DELETE?


def extract_fragments_from_fasta(file, fragments_list):
	with open(file) as fh:
		for sequence in FastaIO.FastaIterator(fh):
			if sequence.id in fragments_list.keys():
				for candidate_peptide in fragments_list[sequence.id]:
					if candidate_peptide[0][0] < candidate_peptide[0][1]:
						fragment = sequence.seq[candidate_peptide[0][0]-1:candidate_peptide[0][1]] # since BLAST is 1-based indexing we need to substract 1
					else:
						fragment = sequence.seq[candidate_peptide[0][1]-1:candidate_peptide[0][0]] # idem
						fragment = fragment.reverse_complement()
					fragments_list[sequence.id][fragments_list[sequence.id].index(candidate_peptide)][0].append(fragment)


def extract_fragments_from_scaffold(scaffold, coordinates):
	if coordinates[0] < coordinates[1]:
		fragment = scaffold[coordinates[0]-1:coordinates[1]]
	else:
		fragment = scaffold[coordinates[1]-1:coordinates[0]]
		fragment = fragment.reverse_complement()
	return(fragment)


def translate_dna(fragments_list, frame = 0):
	for scaffold, scaffold_hits in fragments_list.items():
		for hit in scaffold_hits:
			peptide = hit[0][2][frame:].translate()
			fragments_list[scaffold][fragments_list[scaffold].index(hit)][0].append(peptide)


def align_peptides(protein, fragments_list):
	# aligner = Align.PairwiseAligner(mode='global', substitution_matrix=Align.substitution_matrices.load("BLOSUM62"),
	# 	open_gap_score = -12, extend_gap_score = -5)
	alignments_list = []
	for scaffold, scaffold_hits in fragments_list.items():
		# print(scaffold)
		for hit in scaffold_hits:
			# print(hit)
			for candidate in hit:
				# print(candidate)
				# print(candidate[3])
				# print(peptide)
				with open("seqa.fa", 'w') as fh:
					fh.write('>%s\n%s\n' % (protein.id, protein.seq))
				with open('seqb.fa', 'w') as fh:
					fh.write('>%s.%s.%s\n%s\n' % (scaffold, str(hit[0][0]), str(hit[0][1]), candidate[3]))

				os.system('needle -asequence seqa.fa -bsequence seqb.fa -gapopen 10 -gapextend 1 -outfile pairwise_seqs.fa -aformat fasta')
				alignment = AlignIO.read('pairwise_seqs.fa', 'fasta')
				alignments_list.append(alignment)
	return(alignments_list)

def align_peptides_simple(protein, peptide):
	with open("pairwise_seqs.fa", 'w') as fh:
		fh.write('>%s\n%s\n>fragment_peptide\n%s\n' % (protein.id, protein.seq, peptide))
	os.system('clustalo --infile pairwise_seqs.fa > pairwise_seqs.clustalo.fa')
	alignment = AlignIO.read('pairwise_seqs.clustalo.fa', 'fasta')
	return(alignment)

def get_scaffold_from_fasta(genome, scaffold):
	with open(genome) as fh:
		for sequence in FastaIO.FastaIterator(fh):
			if sequence.id == scaffold:
				return(sequence.seq)

def extend_candidate_peptide(candidate_peptide, scaffold, protein_homolog, direction = 'downstream', order = 0):
	order += 1
	print(candidate_peptide)
	# Get coordinates of next fragment to evaluate (1-based indexing since these are BLAST coordinates)
	if direction == 'downstream':
		if candidate_peptide[0] < candidate_peptide[1]:
			next_fragment_start = candidate_peptide[1] + 1
			next_fragment_end = candidate_peptide[1] + 90
		else:
			next_fragment_start = candidate_peptide[0] + 90
			next_fragment_end = candidate_peptide[0] + 1

	if direction == 'upstream':
		if candidate_peptide[0] < candidate_peptide[1]:
			next_fragment_start = candidate_peptide[0] - 90
			next_fragment_end = candidate_peptide[0] -1
		else:
			next_fragment_start = candidate_peptide[1] - 1
			next_fragment_end = candidate_peptide[1] - 90

	if 0 < next_fragment_start <= len(scaffold) and 0 < next_fragment_end <= len(scaffold) and next:
		next_fragment = extract_fragments_from_scaffold(scaffold, [next_fragment_start, next_fragment_end]) # remember this already does reverse complement if needed
		# print(next_fragment)

		# Translate the fragment to peptide in all 3 reading frames
		peptides_list = [next_fragment.translate(), next_fragment[1:].translate(), next_fragment[2:].translate()]
		# print("Lenght of peptide list is: " + str(len(peptides_list)))

		# Align each peptide to the protein homolog, evaluate the alignment and if it is good, append the fragment
		for index in range(0, len(peptides_list)):
			print('TESTING PEPTIDE:')
			print(peptides_list[index])
			print('PEPTIDE FRAME: ' + str(index))
			alignment = align_peptides_simple(protein = protein_homolog, peptide = peptides_list[index])
			# print(alignment)

			# find where the peptide fragment has the end gaps to not count these in the PID calculation
			within_seq = 0
			identical_positions = 0
			for position in range(0, alignment.get_alignment_length()):
				# print(alignment[0, position])
				# print(alignment[1, position])
				if alignment[1, position] != '-' and not within_seq:
					start = position
					within_seq = 1
				if alignment[1, position] == '-' and alignment[1, position -1] != '-':
					end = position - 1
				if alignment[0, position] == alignment[1, position]:
					identical_positions += 1

			# print(start)
			# print(end)
			# print(identical_positions)

			# Compute PID without considering start/end gaps
			percent_identity = 100*identical_positions/(end - start + 1)
			print('PERCENT IDENTITY:')
			print(percent_identity)

			if percent_identity > 20:
				if index == 0:
					corrected_start = next_fragment_start
					corrected_end = next_fragment_end
				elif index == 1:
					if next_fragment_start < next_fragment_end:
						corrected_start = next_fragment_start + 1
						corrected_end = next_fragment_end - 2
					else:
						corrected_start = next_fragment_start - 1
						corrected_end = next_fragment_end +2
				elif index == 2:
					if next_fragment_start < next_fragment_end:
						corrected_start = next_fragment_start + 2
						corrected_end = next_fragment_end - 1
					else:
						corrected_start = next_fragment_start - 2
						corrected_end = next_fragment_end + 1

				next_fragment_entry = [corrected_start, corrected_end, next_fragment, peptides_list[index]]

				yield (next_fragment_entry, order)
				# yield candidate_peptide
				yield from extend_candidate_peptide(candidate_peptide = next_fragment_entry, scaffold = scaffold, 
					protein_homolog = protein_homolog, direction = direction, order = order)
	# 		else:
	# 			yield candidate_peptide
	# else:
	# 	yield candidate_peptide


if __name__ == '__main__':
	__version__ = "0.0.1"
	
	#### PARSE ARGS ####
	args = docopt(__doc__)
	print(args)

	# if len(args['--genome'].split('/')) > 1:
	# 	args['--genome_dir'] = '/'.join(args['--genome'].split('/')[0:-1])+'/'
	# 	args['--genome'] = args['--genome'].split('/')[-1]

	# if len(args['--proteins'].split('/')) > 1:
	# 	args['--proteins_dir'] = '/'.join(args['--proteins'].split('/')[0:-1])+'/'
	# 	args['--proteins'] = args['--proteins'].split('/')[-1]

	# args['--tblastn_wordsize'] = int(args['--tblastn_wordsize'])
	# args['--tblastn_threads'] = int(args['--tblastn_threads'])
	# args['--tblastn_max_evalue'] = float(args['--tblastn_max_evalue'])

	#### MAIN ####
	if args['--verbose']:
		print("Building BLAST database...")
	
	# Make BLAST nucleotide database with the query genome where we wish to find pseudogenes
	make_blast_db(genome=args['--genome'])

	if args['--verbose']:
		print("Executing tBLASTn search...")

	# Run the tblastn search of the closest homolog proteins against the query genome
	tblastn_output = tblastn(query=args['--proteins'], target=args['--genome'], wordsize=args['--tblastn_wordsize'], 
		matrix=args['--tblastn_matrix'], max_evalue=args['--tblastn_max_evalue'], seg_filter=args['--tblastn_seg_filter'], 
		threads = args['--tblastn_threads'], outprefix=args['--outprefix'])

	# The rest of steps are done protein homolog by protein homolog
	with open(args['--proteins']) as proteins_fh:
		for protein in FastaIO.FastaIterator(proteins_fh):
			protein_id = protein.id
			
			if args['--verbose']:
				print("Processing results for protein %s" % protein_id)

			# Get which positions in the query genome that match the protein homolog
			primary_seeds = find_primary_seed_coordinates(tblastn_output, protein_id)
			
			# print(primary_seeds)

			# Extract the nucleotide sequence corresponding to those positions
			extract_fragments_from_fasta(file=args['--genome'], fragments_list=primary_seeds)

			# print(primary_seeds)

			# Translate the nucleotide sequence to AA
			translate_dna(fragments_list=primary_seeds, frame = 0) # since this comes from tblastn the frist frame is already the good one

			# print(primary_seeds)

			# Align the peptide seed to the protein homolog
			primary_seed_alignment_list = align_peptides(protein = protein, fragments_list=primary_seeds)
			AlignIO.write(primary_seed_alignment_list, 'alignments' + protein.id + '.all.fa', 'fasta')

			# Conduct seed extension
			for scaffold, scaffold_candidate_peptides in primary_seeds.items():
				print('SCAFFOLD: ' + scaffold)
				scaffold_seq = get_scaffold_from_fasta(args['--genome'], scaffold)
				# print(len(scaffold_seq))
				for candidate_peptide in scaffold_candidate_peptides: # here candidate peptides means seeds
					reconstructed_peptides_downstream = []
					reconstructed_peptides_upstream = [[]]
					reconstructed_peptides_downstream.append(candidate_peptide)
					peptide_index = scaffold_candidate_peptides.index(candidate_peptide)
					candidate_peptide
					print('CANDIDATE PEPTIDE SEED: ')
					print(candidate_peptide)
					print('EXTENDING DOWNSTREAM OF THE GENOME:')
					for peptide_tuple in extend_candidate_peptide(candidate_peptide[0], scaffold_seq, protein, direction = 'downstream', order = 0):
						# check the order and based on the length of the already built peptide duplicate or not
						new_peptide = peptide_tuple[0]
						order = peptide_tuple[1]
						if len(reconstructed_peptides_downstream[-1]) == order:
							reconstructed_peptides_downstream[-1].append(new_peptide)
						elif len(reconstructed_peptides_downstream[-1]) == order + 1:
							reconstructed_peptides_downstream.append(reconstructed_peptides_downstream[-1][0:-1] + [new_peptide])

						# print(peptide_tuple)
					print('EXTENDING UPSTREAM OF THE GENOME:')
					for peptide_tuple in extend_candidate_peptide(candidate_peptide[0], scaffold_seq, protein, direction = 'upstream', order = 0):
						new_peptide = peptide_tuple[0]
						order = peptide_tuple[1]
						if len(reconstructed_peptides_upstream[-1]) == order - 1:
							reconstructed_peptides_upstream[-1].insert(0, new_peptide)
						elif len(reconstructed_peptides_upstream[-1]) == order:
							reconstructed_peptides_upstream.append([new_peptide] + reconstructed_peptides_upstream[-1][1:])

					# print(peptide_tuple)
					reconstructed_peptides_complete = []
					for upstream_peptide in reconstructed_peptides_upstream:
						for downstream_peptide in reconstructed_peptides_downstream:
							reconstructed_peptides_complete.append(upstream_peptide + downstream_peptide)

					for index in range(0, len(reconstructed_peptides_complete)):
						this_peptide = reconstructed_peptides_complete[index]
						if this_peptide[0][0] > this_peptide[0][1]:
							reconstructed_peptides_complete[index].reverse()
				
					# Write reconstructed peptides to a file (temporary, output will be formatted as GFF and FASTA later on)
					# print(reconstructed_peptides)
					with open(args['--outprefix'] + '.' + protein.id + '.reconstructed_peptides.txt', 'a') as fh:
						fh.write(scaffold + '\n')
						for peptide in reconstructed_peptides_complete:
							fh.write(str(peptide) + '\n')

	if args['--verbose']:
		print("Execution finished.")

#### END ####