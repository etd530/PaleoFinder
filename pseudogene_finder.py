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
import os
import pandas as pd


#### FUNS ####
def make_blast_db(genome):
	os.system("makeblastdb -dbtype nucl -in " + genome)


def tblastn(query, target, wordsize, matrix, max_evalue, seg_filter, threads, outprefix):
	blast_file = outprefix + ".wordsize" + wordsize + "." + matrix + ".evalue" + max_evalue + ".seg" + seg_filter.replace('"', '') + ".out"
	blast_file = blast_file.replace(' ', '_')
	tblastn_command="tblastn -query " + query + " -db " + target + " -word_size " + wordsize + " -matrix " + matrix + " -evalue " + max_evalue + " -outfmt 6 -num_threads " + threads + " -seg " + seg_filter + " -out " + blast_file
	os.system(tblastn_command)
	tblastn_results = pd.read_csv(blast_file, sep='\t', header = None)
	tblastn_results.rename(columns={0: 'qseqid', 1: 'sseqid', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'}, inplace = True)
	return(tblastn_results)

def find_primary_seed_coordinates(tblastn_df, query_seq_id):
	tblastn_df_subset = tblastn_df[tblastn_df['qseqid'] == query_seq_id]
	primary_seeds = {}
	for index, row in tblastn_df_subset.iterrows():
		newkey = row['sseqid']
		if newkey in primary_seeds.keys():
			primary_seeds[newkey].append((row['sstart'], row['send']))
		else:
			primary_seeds[newkey] = [(row['sstart'], row['send'])]
		
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
			if sequence.id == seqid:
				fragment = sequence.seq[start:end]
	return fragment

# def translate_three_frames(seq, forward = True):
# 	if forward:
		




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
	
	make_blast_db(genome=args['--genome'])

	if args['--verbose']:
		print("Executing tBLASTn search...")

	tblastn_output = tblastn(query=args['--proteins'], target=args['--genome'], wordsize=args['--tblastn_wordsize'], matrix=args['--tblastn_matrix'], 
		max_evalue=args['--tblastn_max_evalue'], seg_filter=args['--tblastn_seg_filter'], threads = args['--tblastn_threads'], outprefix=args['--outprefix'])


	with open(args['--proteins']) as proteins_fh:
		for protein in FastaIO.FastaIterator(proteins_fh):
			protein_id = protein.id
			
			if args['--verbose']:
				print("Processing results for protein %s" % protein_id)

			primary_seeds = find_primary_seed_coordinates(tblastn_output, protein_id)
			
			print(primary_seeds)

			for index in range(0, len(primary_seeds)):
				current_seed = primary_seeds[index]
				current_seed_sequence = extract_fragments_from_fasta(file=args['--genome'], seqid=current_seed[0], start=current_seed[1], end=current_seed[2])
				primary_seeds[index].append(current_seed_sequence)

			

	if args['--verbose']:
		print("Execution finished.")

#### END ####