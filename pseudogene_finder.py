#!/usr/bin/env python3


"""
Reconstruct highly degraded pseudogenes by means of sliding window-based subsequent alignments from a seed alignment.

Usage: 
	pseudogene_finder.py runall --proteins=FASTA --genome=FASTA --blastp_db=STR --parent_taxid=INT [--tblastn_wordsize=INT --tblastn_matrix=STR --tblastn_max_evalue=FLT --tblastn_seg_filter=STR --tblastn_threads=INT --blastp_wordsize=INT --blastp_matrix=STR --blastp_max_evalue=FLT --blastp_threads=INT --outprefix=STR --diamond --diamond_block_size=FLT --excluded_taxids=STR --no_gap_bridging] [-v|--verbose] [-h|--help]
	pseudogene_finder.py tblastn --proteins=FASTA --genome=FASTA [--tblastn_wordsize=INT --tblastn_matrix=STR --tblastn_max_evalue=FLT --tblastn_seg_filter=STR --tblastn_threads=INT --outprefix=STR] [-v|--verbose] [-h|--help]
	pseudogene_finder.py extend --proteins=FASTA --genome=FASTA --tblastn_output=STR [--outprefix=STR --no_gap_bridging] [-v|--verbose] [-h|--help]
	pseudogene_finder.py blastp --blastp_db=STR [--blastp_wordsize=INT --blastp_matrix=STR --blastp_max_evalue=FLT --blastp_threads=INT --diamond --diamond_block_size=FLT] [-v|--verbose] [-h|--help]
	pseudogene_finder.py filter_blastp --proteins=FASTA --blastp_output=STR --parent_taxid=INT [--excluded_taxids=STR --outprefix=STR] [-v|--verbose] [-h|--help]

    Options:
        -p, --proteins FASTA                      Input query proteins for that serve as reference to find pseudogenes, in FASTA format. Can contain multiple sequences
        -g, --genome FASTA                        Input target genome in which pseudogenes are to be found, in FASTA format.
        --tblastn_output STR                      Name of the file containing the tblastn output to use for the extension phase.
        --blastp_db STR                           Path to the blastp/diamond database.
        --tblastn_wordsize INT                    Word size to use for the initial tblastn search [default: 3].
        --tblastn_matrix STR                      Alignment scoring matrix to use for the initial tblastn search [default: BLOSUM62].
        --tblastn_max_evalue FLT                  Maximum e-value to keep a hit in the initial tblastn search [default: 50].
        --tblastn_seg_filter STR                  Parameters for the SEG masking of low complexity regions in the query proteins [default: "10 1.0 1.5"].
        --tblastn_threads INT                     Number of threads to use to run the initial tblastn search [default: 10].
        --no_gap_bridging                         Binary: whether or not to run the gap bridging step.
        --blastp_wordsize INT                     Word size to use for blastp [default: 3].
        --blastp_matrix STR                       Alignment scoring matrix to use for blastp [default: BLOSUM62].
        --blastp_max_evalue FLT                   Maximum e-value to keep a hit in the blastp search [default: 0.2].
        --blastp_threads INT                      Number of threads to use for blastp search [default: 10].
        --blastp_output STR                       Name of the file containing the blastp/diamond output to filter.
        --diamond                                 Use diamond instead of NCBI's blastp for the final validation.
        --diamond_block_size FLT                  Block size to use for diamond; note diamond will use up to about 6 times this value in GBs of RAM [default: 2]
        --parent_taxid INT                        Taxid from NCBI's Taxonomy database specifying a taxa to which the blastp hits are expected to belong. Required to filter the blastp output by taxonomic identity.
        --excluded_taxids STR                     Comma-separated list of taxids to exclude from the blastp results. Usually you want to exlucde the organisms whose genome you are screening to avoid self-hits.
        -o, --outprefix STR                       Prefix to use for output files [default: pseudogene_finder].
        -v, --verbose                             Print the progressions of the program to the terminal (Standard Error).
        -h, --help                                Show this help message and exit.
"""

#### LIBS ####
from docopt import docopt      # to create the argument parser
from Bio.SeqIO import FastaIO  # to work with fasta sequences
from Bio import Align, AlignIO # to work with alignments
import os                      # to send Linux commands
import pandas as pd            # to work with dataframes
import glob                    # to nicely list files from directories
import sys                     # to exit with error messages or fine
import numpy as np             # to do fast array maths
import subprocess              # to send Linux commands and capture output
import taxopy                  # to work with NCBI Taxonomy database
import re                      # to work with regular expressions
import math                    # for mathy stuff. duh

#### FUNS ####
def make_blast_db(genome):
	"""
	Build a BLAST nucleotide database from the input genome.

	Arguments:
		genome: the name of the genome file (in FASTA) from which to build the database.

	Returns:
		Nothing, the database is done by a call to makeblastdb.	
	"""
	if not os.path.exists(genome+'.nhr'):
		os.system("makeblastdb -dbtype nucl -in " + genome)
	else:
		print("WARNING: BLAST database already found for " + genome + ". Skipping database building.") # Could implement a force flag if they want to rerun, but let's see if we end up using submodules...

def tblastn(query, target, wordsize, matrix, max_evalue, seg_filter, threads, outprefix):
	"""
	Do a tblastn search of the selected homologous proteins against the database built from the target genome.

	Arguments:
		query: path of the file containing the query proteins for the tblastn.
		target: path of the genome from which the database was built.
		wordsize: word size to use for the tblastn search.
		matrix: aminoacid substitution matrix to use to score the local alignments.
		max_evalue: maximum evalue to report a hit.
		seq_filter: parameters for the SEG algorithm to mask low complexity regions in the query proteins.
		threads: number of threads to use for the tblastn search.
		outprefix: prefix to use to build the name for the tblastn output file.

	Returns:
		A pandas dataframe with the tblastn results.
	"""
	blast_file = outprefix + ".wordsize" + wordsize + "." + matrix + ".evalue" + max_evalue + ".seg" + seg_filter.replace('"', '') + ".out"
	blast_file = blast_file.replace(' ', '_')

	if not os.path.exists(blast_file):
		tblastn_command="tblastn -query " + query + " -db " + target + " -word_size " + wordsize + " -matrix " + matrix + " -evalue " + max_evalue + " -outfmt 6 -num_threads " + threads + " -seg " + seg_filter + " -out " + blast_file
		os.system(tblastn_command)
	else:
		print("WARNING: tBLASTn results found at " + blast_file + ". Skipping tBLASTN step.")
	
	try:
		tblastn_results = pd.read_csv(blast_file, sep='\t', header = None)
	except pd.errors.EmptyDataError:
		sys.exit("No hits obtained from tblastn, exiting program.")

	tblastn_results.rename(columns={0: 'qseqid', 1: 'sseqid', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'}, inplace = True)
	
	return(tblastn_results)

def find_primary_seed_coordinates(tblastn_df, query_seq_id):
	"""
	Find the coordinates of the primary seeds (i.e. blast hits) form the blast output.

	Arguments:
		tblastn_df: Pandas dataframe containing the tblastn results.
		query_seq_id: ID of the query protein for which to recover the hit coordinates in the target genome.

	Returns:
		A dictionary containing the name of the scaffolds of the genome as keys, and a list with the start and end positions for each of the hits found in that scaffold.
	"""
	tblastn_df_subset = tblastn_df[tblastn_df['qseqid'] == query_seq_id]
	primary_seeds = {}
	for index, row in tblastn_df_subset.iterrows():
		newkey = row['sseqid']
		if newkey in primary_seeds.keys():
			primary_seeds[newkey].append([[row['sstart'], row['send'], row['qstart'], row['qend']]])
		else:
			primary_seeds[newkey] = [[[row['sstart'], row['send'], row['qstart'], row['qend']]]]
		
		# print(newkey)
		# primary_seeds

	# = list(tblastn_df_subset.apply(lambda row: [row['sseqid'], row['sstart'], row['send']], axis=1))
	return(primary_seeds)

def extract_fragments_from_fasta(file, fragments_list):
	"""
	Extract the sequences of the tblastn hits from the genome sequence.

	Arguments:
		file: the FASTA file containing the target genome.
		fragments_list: a dictionary containing the genome scaffolds' names as keys and a list with the coordinates of the hits in that scaffold as value.

	Returns:
		Nothing, it appends the sequences to the list.
	"""
	with open(file) as fh:
		for sequence in FastaIO.FastaIterator(fh):
			if sequence.id in fragments_list.keys():
				for candidate_peptide in fragments_list[sequence.id]:
					if candidate_peptide[0][0] < candidate_peptide[0][1]:
						fragment = sequence.seq[candidate_peptide[0][0]-1:candidate_peptide[0][1]] # since BLAST is 1-based indexing we need to substract 1
					else:
						fragment = sequence.seq[candidate_peptide[0][1]-1:candidate_peptide[0][0]] # idem
						fragment = fragment.reverse_complement()
					fragments_list[sequence.id][fragments_list[sequence.id].index(candidate_peptide)][0].insert(2, fragment)

def extract_fragments_from_scaffold(scaffold, coordinates):
	"""
	Extract a sequence from a scaffold given as sequence object.

	Arguments:
		scaffold: the sequence object from which to recover the subsequences.
		coordinates: the pair of coordinates of a given hit.

	Returns:
		The sequence fragment, in the strand that gets translated.
	"""
	if coordinates[0] < coordinates[1]:
		# print('Peptide comes from the forward strand so no need to reverse complement.')
		fragment = scaffold[coordinates[0]-1:coordinates[1]]
	else:
		# print('Peptide comes from the reverse strand so oligonucleotide will be reverse complemented.')
		fragment = scaffold[coordinates[1]-1:coordinates[0]]
		fragment = fragment.reverse_complement()
	return(fragment)

def translate_dna(fragments_list, frame = 0):
	"""
	Translate a DNA sequence to AA sequence in the specified reading frame.

	Arguments:
		fragments_list: a dictionary with scaffold names as keys and a list containing the sequence to translate as the third element as values.
		frame: the desired reading frame to consider for translation.

	Returns:
		Nothing, the AA sequence is appended into the list.
	"""
	for scaffold, scaffold_hits in fragments_list.items():
		for hit in scaffold_hits:
			peptide = hit[0][2][frame:].translate()
			fragments_list[scaffold][fragments_list[scaffold].index(hit)][0].insert(3, peptide)

def align_peptides_simple(protein, peptide):
	"""
	Align 2 aminoacid sequences using the specified program (either EMBOSS Needle for global alignment of lalign36 for non-overlapping local alignments).

	Arguments:
		protein: one of the AA sequences.
		peptide: the other AA sequence.
	Returns:
		An alignment object.
	"""
	# Check peptide size, if it is 30AA do global alignment, else it should be 100AA, do local alignment
	# print(peptide)
	if len(peptide) == 30 or len(peptide) == 29: #putting just if len() < 31 should do, but I want to make sure there is no unforeseen stuff; change later?
		# Align with EMBOSS-Needle
		with open("seqa.fa", 'w') as fh:
			fh.write('>%s\n%s\n' % (protein.id, protein.seq))
		with open("seqb.fa", "w") as fh:
			fh.write('>fragment_peptide\n%s\n' % peptide)
		os.system('needle -asequence seqa.fa -bsequence seqb.fa -gapopen 10 -gapextend 1 -outfile pairwise_seqs.fa -aformat fasta -auto Y -sprotein1 Y -sprotein2 Y')
		status = os.system('cat pairwise_seqs.fa >> alignments/pairwise_seqs.tmp.fa')
		alignment = AlignIO.read('pairwise_seqs.fa', 'fasta')
		# Align with Clustal Omega (Needle works best for our case for now)
		# with open("input_seqs.fa", 'w') as fh:
		# 	fh.write('>%s\n%s\n>fragment_peptide\n%s\n' % (protein.id, protein.seq, peptide))
		# os.system('clustalo --infile input_seqs.fa > pairwise_seqs.fa')
		# os.system('cat pairwise_seqs.fa >> alignments/pairwise_seqs.tmp.fa')
		# alignment = AlignIO.read('pairwise_seqs.fa', 'fasta')
	elif len(peptide) == 100 or len(peptide) == 99: #same as above, putting just len()<101 should do
		# Align with lalign36
		with open("seqa.fa", 'w') as fh:
			fh.write('>%s\n%s\n' % (protein.id, protein.seq))
		with open("seqb.fa", "w") as fh:
			fh.write('>fragment_peptide\n%s\n' % peptide)
		# NOTE: THIS LINE WILL NOT WORK WITHOUT USERS ADDING THE PROGRAM TO THE PATH, NEED TO SEE HOW TO FIX THIS
		os.system('lalign36 -O lalign.aln -m 3 -3 -C 20 -Q -q seqa.fa seqb.fa > lalign.log')
		status = os.system('num1=`egrep -c "^>" lalign.aln`; num2=2; [ $num1 -eq $num2 ]')
		if not status:
			os.system('''cat lalign.aln | head -n -10 | tail -n +20 > pairwise_seqs.fa && \
				seqindex1=`grep -n '^>' pairwise_seqs.fa | cut -f1 -d':' | head -n1` && \
				seqindex2=`grep -n '^>' pairwise_seqs.fa | cut -f1 -d':' | tail -n1` && seqname=`head -n1 seqa.fa | sed -E 's/\//_/g'` && \
				sed -Ei "${seqindex1}s/.*/${seqname}/" pairwise_seqs.fa && seqname=`head -n1 seqb.fa | sed -E 's/\//_/g'` && \
				sed -Ei "${seqindex2}s/.*/${seqname}/" pairwise_seqs.fa && cat pairwise_seqs.fa >> alignments/pairwise_seqs.tmp.fa''')
	# if os.path.isfile('pairwise_seqs.fa'):
	if not status:
		alignment = AlignIO.read('pairwise_seqs.fa', 'fasta')
		return(alignment)

def align_bridges(protein_gap, bridge):
	"""
	Align 2 aminoacid sequences using EMBOSS Needle.

	Arguments:
		protein_gap: the gap in the protein homolog that needs to be covered.
		bridge: the peptide obtained by translating the genome that should be able to cover the gap in the protein homolog.
	Returns:
		An alignment object.
	"""
	# Align with EMBOSS-Needle
	with open("seqa.fa", 'w') as fh:
		fh.write('>%s\n%s\n' % (protein_gap.id, protein_gap.seq))
	with open("seqb.fa", "w") as fh:
		fh.write('>bridge_peptide\n%s\n' % bridge)
	os.system('needle -asequence seqa.fa -bsequence seqb.fa -gapopen 10 -gapextend 1 -outfile pairwise_seqs.fa -aformat fasta -auto Y -sprotein1 Y -sprotein2 Y')
	status = os.system('cat pairwise_seqs.fa >> alignments/pairwise_seqs.tmp.fa')
	alignment = AlignIO.read('pairwise_seqs.fa', 'fasta')
	return(alignment)

def get_scaffold_from_fasta(genome, scaffold):
	"""
	Obtain the sequence of a given scaffold from a genome FASTA.

	Arguments:
		genome: the FASTA file containing the genome.
		scaffold: the name of the scaffold whose sequence is to be recovered.

	Returns:
		The sequence of the scaffold.
	"""
	with open(genome) as fh:
		for sequence in FastaIO.FastaIterator(fh):
			if sequence.id == scaffold:
				return(sequence.seq)

def extend_candidate_peptide(candidate_peptide, scaffold, protein_homolog, direction = 'downstream', order = 0, fragment_size = 90, do_local = True, reading_frame = 'All'):
	"""
	Given a short AA sequence, extend it based on similarity of subsequent fragments from the genome of origin to a reference protein sequence.

	Arguments:
		candidate_peptide: a list containing the start and end coordinates of the peptide in the genome's scaffold and the nucleotide sequence.
		scaffold: a Seq object containing the sequence of the scaffold form which the peptide was translated.
		protein_homolog: the protein to which compare subsequent fragments of the protein.
		direction: string indicating whether the peptide must be extended upstream or downstream on the scaffold.
		order: number indicating the number of iterations, used internally to recover the order of the peptides.
		fragment_size: the size of the fragments to use for the successive alignments to extend the peptide.
		do_local: whether a local alignment of a longer fragment should be attempted when all 90bp fragments align poorly.
	Yields:
		a AA fragment each time that one passes the established thresholds of similarity and contiguity.
	"""

	order += 1
	print('INITIAL GENOMIC COORDINATES:')
	print(candidate_peptide[0])
	print(candidate_peptide[1])
	print('FRAGMENT SIZE:')
	print(fragment_size)
	# Get coordinates of next fragment to evaluate (1-based indexing since these are BLAST coordinates)
	if direction == 'downstream':
		if candidate_peptide[0] < candidate_peptide[1]:
			next_fragment_start = candidate_peptide[1] + 1
			next_fragment_end = candidate_peptide[1] + fragment_size
		else:
			next_fragment_start = candidate_peptide[0] + fragment_size
			next_fragment_end = candidate_peptide[0] + 1

	if direction == 'upstream':
		if candidate_peptide[0] < candidate_peptide[1]:
			next_fragment_start = candidate_peptide[0] - fragment_size
			next_fragment_end = candidate_peptide[0] -1
		else:
			next_fragment_start = candidate_peptide[1] - 1
			next_fragment_end = candidate_peptide[1] - fragment_size
	print('INITIAL COORDINATES OF NEXT FRAGMENT:')
	print(next_fragment_start)
	print(next_fragment_end)
	# check if the next fragment would be within the scaffold, and not outside (i.e. if the scaffold has at least 90/300bp left)
	# UPDATE THIS TO CREATE A SHORTER FRAGMENT INSTEAD OF JUST NOT DOING ANYTHING WIHT THE SCAFFOLD EXTREMES
	if 0 < next_fragment_start <= len(scaffold) and 0 < next_fragment_end <= len(scaffold) and next:
		next_fragment = extract_fragments_from_scaffold(scaffold, [next_fragment_start, next_fragment_end]) # remember this already does reverse complement if needed
		print('NEXT FRAGMENT:')
		print(next_fragment)

		# Translate the fragment to peptide in all 3 reading frames, except for the local alignment which needs to be done only for the same frame as the one that failed with global alignment
		if reading_frame == 'All':
			peptides_list = [next_fragment.translate(), next_fragment[1:].translate(), next_fragment[2:].translate()]
			print("Length of peptide list is: " + str(len(peptides_list)))
		else:
			peptides_list = [next_fragment[reading_frame:].translate()]

		# Align each peptide to the protein homolog, evaluate the alignment and if it is good, append the fragment
		for index in range(0, len(peptides_list)):
			print('TESTING PEPTIDE:')
			print(peptides_list[index])
			print('PEPTIDE LENGTH:')
			print(len(peptides_list[index]))
			if len(peptides_list) > 1:
				print('PEPTIDE FRAME: ' + str(index))
			else:
				print('PEPTIDE FRAME: ' + str(reading_frame))
			print('PEPTIDE ORDER: ' + str(order))
			alignment = align_peptides_simple(protein = protein_homolog, peptide = peptides_list[index])
			# print(alignment)

			# make sure you actually got an alignment before continuing (local alignment may not return anything)
			if alignment is not None:
				# find where the peptide fragment has the end gaps to not count these in the PID calculation
				# also find the start and end coordinates relative to the homologous protein
				# FOR NOW THIS DOES NOT ACCOUNT FOR GAPS IN THE ALIGNMENT IN THE PROTEIN HOMOLOG
				within_seq = 0
				identical_positions = 0
				position_in_homolog = 0
				position_in_fragment = 0
				end_found = False
				print('ALIGNMENT:')
				print(alignment)
				for position in range(0, alignment.get_alignment_length()):
					if alignment[0, position] != '-':
						position_in_homolog += 1
						if alignment[1, position] != '-':
							position_in_fragment += 1
							end_found = False # set end found to False each time you find an AA, to prevent internal gaps from being considered the end if there is not end gap
							if not within_seq:
								start = position
								homolog_start = position_in_homolog
								fragment_start = position_in_fragment - 1
								within_seq = 1
						if within_seq and alignment[1, position] == '-' and alignment[1, position -1] != '-': # this only founds the end if there is an end gap!
							# print("END FOUND")
							# print("PREVIOUS POSITION:")
							# print(alignment[1, position -1 ])
							end = position - 1
							end_found = True
							homolog_end = position_in_homolog - 1
							fragment_end = position_in_fragment -1
						if alignment[0, position] == alignment[1, position]:
							identical_positions += 1
					elif alignment[1, position] != '-': # if there is gap in homolog, still check also if there is gap in fragment to add to the position if needed
						position_in_fragment += 1
				# after checking all positions, first make sure that start and homolog_start have been found, else it means there is not overlap in the seqs (i.e. no real alignment, all is gaps)
				if 'homolog_start' in locals() and 'start' in locals():
					# if end has not been found because there is no end gap, it means the end is the end of the actual alignment
					if not end_found:
						end = alignment.get_alignment_length() - 1
						homolog_end = position_in_homolog - 1
						fragment_end = position_in_fragment - 1

					# Check if candidate fragment aligns fully (it is "bounded" by the reference homolog); if not, mark as outside the bounds
					print("ORIGINAL next fragment:")
					print(peptides_list[index])
					print("COORDINATES OF ALIGNED PART IN FRAGMENT: %s - %s" % (str(fragment_start), str(fragment_end)))
					if len(peptides_list[index]) > (end - start + 1):
						print("A part of the fragment goes over the edges of the reference homolog, will keep only the part that aligns.")
						outbound = True
					else:
						print("Peptide candidate is fully bounded by the reference homolog.")
						outbound = False
					# print("SHORTENED next fragment:")
					# print(next_fragment)

					# if doing local alignment, the output only contains the subsequences that align correctly, so we need to locate the homolog subsequence in the sequence of the full homolog and transform the coordinates
					if fragment_size == 300:
							aligned_region = alignment[0].seq
							# print(alignment[1].seq.replace('-', ''))
							homolog_start = protein_homolog.seq.find(aligned_region.replace('-', '')) + 1 # add one since we are doing 1-based indexing
							homolog_end = homolog_start + len(aligned_region) - 1

					# check if current fragment aligned is contiguous with the previously aligned one
					# print(candidate_peptide[4])
					# print(candidate_peptide[5])
					# print(homolog_start)
					# print(homolog_end)
					if direction == 'downstream': # if fragments go downstream of the genome lead strand
						if candidate_peptide[0] < candidate_peptide[1]: # if the gene is in the lead strand
							if 0 < homolog_start - candidate_peptide[5] <= 10:
								contiguous = True
							else:
								contiguous = False
						else:
							if 0 < candidate_peptide[4] - homolog_end <= 10:
								contiguous = True
							else:
								contiguous = False
					else:
						if candidate_peptide[0] < candidate_peptide[1]:
							if 0 < candidate_peptide[4] - homolog_end <= 10:
								contiguous = True
							else:
								contiguous = False
						else:
							if 0 < homolog_start - candidate_peptide[5] <= 10:
								contiguous = True
							else:
								contiguous = False

					# Compute PID without considering start/end gaps
					# print(start)
					# print(end)
					percent_identity = 100*identical_positions/(end - start + 1)
					
					print('PERCENT IDENTITY:')
					print(percent_identity)
					print('START IN FRAGMENT:')
					print(fragment_start)
					print('END IN FRAGMENT:')
					print(fragment_end)

					# If the fragment passes the thresholds, get the final coordinates and create the fragment entry for the list
					if contiguous and percent_identity > 20:
						# First we correct the coordinates for the reading frame
						# *IMPORTANT* NOTE: we need the OR condition because index is always zero when trying local alignment since we only translate in one frame in that case!!!
						if reading_frame == 'All' and index == 0 or reading_frame == 0:
							print('READING FRAME 0')
							corrected_start = next_fragment_start
							corrected_end = next_fragment_end
						elif reading_frame == 'All' and index == 1 or reading_frame == 1:
							print('READING FRAME 1')
							if next_fragment_start < next_fragment_end:
								corrected_start = next_fragment_start + 1
								corrected_end = next_fragment_end - 2
							else:
								corrected_start = next_fragment_start - 1
								corrected_end = next_fragment_end +2
						elif reading_frame == 'All' and index == 2 or reading_frame == 2:
							print('READING FRAME 2')
							if next_fragment_start < next_fragment_end:
								corrected_start = next_fragment_start + 2
								corrected_end = next_fragment_end - 1
							else:
								corrected_start = next_fragment_start - 2
								corrected_end = next_fragment_end + 1
						# print('CODING FRAME:')
						# print(index)
						# print('CORRECTED GENOMIC COORDINATES FOR THE FRAME:')
						# print(corrected_start)
						# print(corrected_end)
						# if fragment is outside bounds of reference homolog, we need to keep ony the part that aligns and its coordinates
						# then if fragment size is 300, it means we did a local alignment, so we need the coordinates for the part that got aligned in the homolog, as well as extracting the peptide subsequence
						# and if fragment is outside the bounds of the reference homolog, we need also to get the coordinates of the aligned framgnet, but it works a bit differently
						
						if fragment_size == 300 or outbound:
							# extract peptide subsequence and find position of aligned region relative to the full fragment
							if fragment_size == 300:
								print('FRAGMENT WAS ASSESSED BY LOCAL ALIGNMENT, EXTRACTING ALIGNED REGION')
								aligned_region = alignment[1].seq
								aligned_region = aligned_region.replace('-', '')
								print('ALIGNED REGION:')
								print(aligned_region)
								aligned_region_start = peptides_list[index].find(aligned_region)
								aligned_region_end = aligned_region_start + len(aligned_region) - 1
							else:
								print('FRAGMENT FROM GLOBAL ALIGNMENT BUT OUTSIDE BOUNDS OF HOMOLOG, EXTRACTING ALIGNED REGION')
								aligned_region = peptides_list[index][fragment_start:fragment_end + 1]
								print('ALIGNED REGION:')
								print(aligned_region)
								aligned_region_start = fragment_start
								aligned_region_end = fragment_end
								print(aligned_region_start)
								print(aligned_region_end)
								print('COORDINATES OF ALIGNED REGION IN ORIGINAL FRAGMENT:')
								print(fragment_start)
								print(fragment_end)

							peptide_head_gap = aligned_region_start
							aligned_region_len = len(aligned_region)
							peptide_end_gap = len(peptides_list[index]) - peptide_head_gap - aligned_region_len
							aligned_nucleotide_len = len(aligned_region)*3
							print('LENGTH OF PORTION OF PEPTIDE THAT ALIGNS:')
							print(aligned_region_len)
							print('LENGTH OF AMINOACIDS NOT ALIGNED IN THE N-TERMINAL:')
							print(peptide_head_gap)
							print('LENGTH OF AMINOACIDS NOT ALIGNED IN THE C-TERMINAL:')
							print(peptide_end_gap)

							# Correct again the genomic coordinates considering only the short fragment that aligns of the initial 300bp; needs to accound for sense/antisense
							if candidate_peptide[0] < candidate_peptide[1]:
								print("PEPTIDE IS ON THE SENSE STRAND:")
								corrected_fragment_len = corrected_end - corrected_start + 1
								nucleotide_head_gap = 3*peptide_head_gap
								print(nucleotide_head_gap)
								corrected_start = corrected_start + nucleotide_head_gap
								print('CORRECTED START IN GENOME:')
								print(corrected_start)
								corrected_end = corrected_start + aligned_nucleotide_len - 1
								print('CORRECTED END IN GENOME:')
								print(corrected_end)
								aligned_nucleotide_start = corrected_start - next_fragment_start # note that we should add +1 BUT since we use this to subset for Python, and which is 0-based and we are using 1-based, we would have to subtract 1, so we save that step
								print('CORRECTED START IN FRAGMENT:')
								print(aligned_nucleotide_start)
								aligned_nucleotide_end = aligned_nucleotide_start + aligned_nucleotide_len - 1 # same as previous line
								print('CORRECTED END IN FRAGMENT:')
								print(aligned_nucleotide_end)
								next_fragment_aligned_region = next_fragment[aligned_nucleotide_start:aligned_nucleotide_end + 1]
								print('ALIGNED REGION FROM GENOME:')
								print(next_fragment_aligned_region)

							else:
								print("PEPTIDE IS ON THE ANTI-SENSE STRAND:")
								nucleotide_head_gap = 3*peptide_end_gap
								corrected_end = corrected_end + nucleotide_head_gap
								corrected_start = corrected_end + aligned_nucleotide_len - 1
								print("CORRECTED START IN GENOME")
								print(corrected_start)
								print("CORRECTED END IN GENOME")
								print(corrected_end)
								print("PREVIOUS END IN GENOME")
								print(next_fragment_end)
								aligned_nucleotide_start = next_fragment_start - corrected_start # note that we should add +1 BUT since we use this to subset for Python, and which is 0-based and we are using 1-based, we would have to subtract 1, so we save that step
								# aligned_nucleotide_end = corrected_end - next_fragment_end # same as previous line
								aligned_nucleotide_end = aligned_nucleotide_start + aligned_nucleotide_len -1
								print("CORRECTED START IN FRAGMENT")
								print(aligned_nucleotide_start)
								print("CORRECTED END IN FRAGMENT")
								print(aligned_nucleotide_end)
								next_fragment_aligned_region = next_fragment[aligned_nucleotide_start:aligned_nucleotide_end + 1]
							print("LENGTH OF NUCLEOTIDES NOT ALIGNED IN THE 5':")
							print(nucleotide_head_gap)
							print("GENOMIC COORDINATES CORRECTED FOR THE ALIGNED FRAGMENT:")
							print(corrected_start)
							print(corrected_end)
							print("DNA FRAGMENT CORRESPONDING TO THE ALIGNED PORTION OF THE PEPTIDE:")
							print(next_fragment_aligned_region)
							print("TRANSLATION TO CHECK IT WORKS OK:")
							print(next_fragment_aligned_region.translate())
						
						# then build the list with the attributes of the fragment to yield # NEED TO UPDATE THIS TO GIVE START AND END INVERTED IF ON ANTISENSE STRAND
						if fragment_size == 90 and not outbound:
							print('FRAGMENT FROM GLOBAL ALIGNMENT AND INSIDE BOUNDS')
							aligned_nucleotide_len = 3*len(peptides_list[index])
							if next_fragment_start < next_fragment_end:
								aligned_nucleotide_start = corrected_start - next_fragment_start
							else:
								aligned_nucleotide_start = next_fragment_start - corrected_start
							
							aligned_nucleotide_end = aligned_nucleotide_start + aligned_nucleotide_len - 1
							next_fragment_aligned_region = next_fragment[aligned_nucleotide_start:aligned_nucleotide_end + 1]
							print(corrected_start)
							print(corrected_end)
							print(next_fragment_start)
							print(next_fragment_end)
							print(next_fragment_aligned_region)
							print(peptides_list[index])
							
							# reverse complement if needed (since we are storing the sequence in the forward strand always), then store fragment
							if next_fragment_start < next_fragment_end:
								print("Building new fragment on a forwrad strand")
								next_fragment_entry = [corrected_start, corrected_end, next_fragment_aligned_region, peptides_list[index], homolog_start, homolog_end]
							else:
								print("Building new fragment on a reverse strand")
								next_fragment_entry = [corrected_start, corrected_end, next_fragment_aligned_region.reverse_complement(), peptides_list[index], homolog_start, homolog_end]
						else:
							if next_fragment_start < next_fragment_end:
								print("Building new fragment on a forwrad strand from local alignment (or outbound global)")
								next_fragment_entry = [corrected_start, corrected_end, next_fragment_aligned_region, aligned_region, homolog_start, homolog_end]
							else:
								print("Building new fragment on a reverse strand form local alignment (or outbound global)")
								next_fragment_entry = [corrected_start, corrected_end, next_fragment_aligned_region.reverse_complement(), aligned_region, homolog_start, homolog_end]

						# check to make sure nucleotide translates correctly to peptide, then yield
						if next_fragment_entry[0] < next_fragment_entry[1]:
							print("Final check in forward strand case.")
							print(next_fragment_entry[2])
							print(next_fragment_entry[2].translate())
							print(next_fragment_entry[3])
							assert next_fragment_entry[2].translate() == next_fragment_entry[3]
						else:
							print("Final check in reverse strand case.")
							print(next_fragment_entry[2])
							print(next_fragment_entry[2].reverse_complement().translate())
							print(next_fragment_entry[3])
							assert next_fragment_entry[2].reverse_complement().translate() == next_fragment_entry[3]
						yield (next_fragment_entry, order)
						
						# lastly, call the next iteration 
						yield from extend_candidate_peptide(candidate_peptide = next_fragment_entry, scaffold = scaffold, 
							protein_homolog = protein_homolog, direction = direction, order = order, do_local = do_local)

			 		# If it does not pass the thresholds AND we have NOT done local alignment yet, start the 300bp local alignment step
					elif do_local:
						print("Peptide not good. Starting round of local alignment")
						yield from extend_candidate_peptide(candidate_peptide = candidate_peptide, scaffold = scaffold,
							protein_homolog = protein_homolog, direction = direction, order = order - 1, fragment_size = 300, do_local = False, reading_frame = index)
					else:
						print("Peptide not good and local alignment tried. Going back one order to test remaining peptides of previous order.")
						pass
				
				# if homolog start is not assigned, it means no positions were aligned at all, do 300 step IF not tried yet
				elif do_local:
					print("Peptide not good. Starting round of local alignment")
					yield from extend_candidate_peptide(candidate_peptide = candidate_peptide, scaffold = scaffold,
						protein_homolog = protein_homolog, direction = direction, order = order - 1, fragment_size = 300, do_local = False, reading_frame = index)
				else:
					print("Peptide not good and local alignment tried. Going back one order to test remaining peptides of previous order.")
					pass

			else:
				print("No local alignment produced for this sequence, it will not be added to the peptide")
				pass
	else:
		# print('INFO: Skipping this peptide, it would be outside the bounds of the scaffold. If this happened with the seed, there will be no alignment output.')
		pass

def number_of_substitutions(seq1, seq2, no_gaps = False):
	"""
	Compute the number of differences between 2 input sequences.

	Arguments:
		seq1: the first sequence.
		seq2: the second sequence.
		no_gaps: whether to consider gaps (indels) as a difference or not. If false, gap positions are ignored

	Returns:
		An integer indicating the number of differences between both sequences.
	"""
	assert len(seq1) == len(seq2)
	if no_gaps:
		return sum(1 if nuc1 != nuc2 and nuc1 != '-' and nuc2 != '-' else 0 for(nuc1, nuc2) in zip(seq1, seq2))
	else:
		return sum(0 if nuc1 == nuc2 else 1 for (nuc1, nuc2) in zip(seq1, seq2))

def gap_bridging(reconstructed_peptides, scaffold, homolog):
	"""
	Given a list of reconstructed candidate peptides, for those with a gap between the fragments both in the genome and when aligned to the protein homolog, try to add more nucleotides of this gap that fill the protein gap.

	Arguments:
		reconstructed_peptides: list containing the reconstructed peptides from the extension step.
		scaffold: the scaffold in which the current candidate peptides are located in the query genome.
		homolog: Seq object of the homolog protein which serves as template to find the candidate pseudogenes

	Returns:
		The list of reconstructed peptides with the gaps reduced/fragments bridged.
	"""
	def build_bridged_fragment(current_fragment, scaffold, bridge_coordinates, forward, homolog_gap):
		"""
		Extend a sequence fragment based on coordinates to extend it.

		Arguments:
		current_fragment: the fragment to extend.
		scaffold: the scaffold from which to get the sequence to extend the fragment.
		bridge_coordinates: the coordinates in the scafofld from which to get the sequence.
		forwards: whether the sequence is on the forward or reverse strand.
		homolog_gap: Seq object containing the sequence gap that needs to be covered in the protein homolog.

		Returns:
			The extended fragment.
		"""
		if forward:
			dna_bridge = extract_fragments_from_scaffold(scaffold, bridge_coordinates)
			# print(current_fragment[2])
			extended_nt_sequence = dna_bridge + current_fragment[2]
			aa_bridge = dna_bridge.translate()
			print("AA BRIDGE:")
			print(aa_bridge)
			print("HOMOLOG GAP:")
			print(homolog_gap)
			assert len(aa_bridge) == homolog_gap
			extended_aa_sequence = aa_bridge + current_fragment[3]
			current_fragment[3] = extended_aa_sequence
			current_fragment[2] = extended_nt_sequence
			current_fragment[0] = bridge_coordinates[0]
			current_fragment[4] = current_fragment[4] - homolog_gap

		else:
			dna_bridge = extract_fragments_from_scaffold(scaffold, bridge_coordinates)
			extended_nt_sequence = current_fragment[2].seq() + dna_bridge
			aa_bridge = dna_bridge.translate()
			assert len(aa_bridge) == homolog_gap
			extended_aa_sequence = current_fragment[3].seq() + aa_bridge
			current_fragment[3] = extended_aa_sequence
			current_fragment[2] = extended_nt_sequence
			current_fragment[1] = bridge_coordinates[1]
			current_fragment[5] = current_fragment[5] + homolog_gap
		return current_fragment

	print('STARTING GAP BRIDGING:')
	bridge_coordinates_list = []
	for peptide_index in range(0, len(reconstructed_peptides)):
		candidate_peptide = reconstructed_peptides[peptide_index]
		print('CURRENT PEPTIDE:')
		print(candidate_peptide)
		for fragment_index in range(0, len(candidate_peptide)):
			current_fragment = candidate_peptide[fragment_index]
			print("CURRENT FRAGMENT:")
			print(current_fragment)
			# for each fragment we look at the gap before it so for the first one there is nothing to do except grab its end positions
			if fragment_index != 0:
				current_start_in_homolog = current_fragment[4]
				if forward:
					current_start_in_genome = current_fragment[0]
				else:
					current_start_in_genome = current_fragment[1]
				print('Start of the fragment in genome: %s' % str(current_start_in_genome))
				print('End of previous fragment in genome: %s' % str(previous_end_in_genome))
				genome_gap = current_start_in_genome - previous_end_in_genome - 1
				homolog_gap = current_start_in_homolog - previous_end_in_homolog -1
				if genome_gap>=3 and homolog_gap:
					print('Both genome and homolog have a gap, with genome gap of at least 3AA.')
					homolog_gap_seq = extract_fragments_from_scaffold(homolog, [previous_end_in_homolog +1, current_start_in_homolog -1])
					if 3*homolog_gap == genome_gap:
						print('DNA gap size matches that of aminoacid gap times 3. Trying to fit the entire fragment in the gap.')
						bridge_start = previous_end_in_genome + 1
						bridge_end = current_start_in_genome -1
						reconstructed_peptides[peptide_index][fragment_index] = build_bridged_fragment(current_fragment, 
							scaffold, [bridge_start, bridge_end], forward, homolog_gap)

					elif genome_gap > 3*homolog_gap:
						print('DNA gap is larger than aminoacid gap times 3, trying to fit a sliding window of the DNA into the homolog peptide.')
						current_candidate_fragment = None
						for i in range(0, genome_gap):
							sub_bridge_start = previous_end_in_genome + i # note this works bc we already took the start and end relative to the fowrard strand
							sub_bridge_end = sub_bridge_start + 3*homolog_gap - 1
							# print(sub_bridge_start)
							# print(sub_bridge_end)
							# print(homolog_gap)
							assert (sub_bridge_end - sub_bridge_start + 1) == 3*homolog_gap
							if sub_bridge_end <= current_start_in_genome - 1:
								new_candidate_fragment = extract_fragments_from_scaffold(scaffold, [sub_bridge_start, sub_bridge_end])
								new_candidate_fragment = [sub_bridge_start, sub_bridge_end, new_candidate_fragment, new_candidate_fragment.translate(), previous_end_in_homolog + 1, current_start_in_homolog - 1]
								# print("NEW CANDIDATE FRAGMENT:")
								# print(new_candidate_fragment)
								# build_bridged_fragment(current_fragment, scaffold, [sub_bridge_start, sub_bridge_end], forward, homolog_gap_seq)
								if current_candidate_fragment is None:
									current_candidate_fragment = new_candidate_fragment
									print(homolog_gap_seq)
									print(current_candidate_fragment)
									subs_num = number_of_substitutions(homolog_gap_seq, current_candidate_fragment[3])
									current_candidate_pid = subs_num/homolog_gap
								else:
									subs_num = number_of_substitutions(homolog_gap_seq, new_candidate_fragment[3])
									new_candidate_pid = subs_num/homolog_gap
									if new_candidate_pid > current_candidate_pid:
										current_candidate_fragment = new_candidate_fragment
										current_candidate_pid = new_candidate_pid
						print('FINISHED CHECKING SLIDING WINDOW. FINAL SUBBRIDGE IS:')
						print(current_candidate_fragment)
						print('RECONSTRUCTED PEPTIDES BEFORE:')
						print(reconstructed_peptides)
						reconstructed_peptides[peptide_index].insert(fragment_index, current_candidate_fragment)
						print('RECONSTRUCTED PEPTIDES AFTER:')
						print(reconstructed_peptides)
						
					else:
						print('DNA gap is smaller than aminoacid gap times 3, trying to cut from the previous fragments and adjust the reading frames to fit a bridge.')
						gaps_difference = 3*homolog_gap - genome_gap
						# Get the bridge to fill the protein gap
						print("Cutting from previous fragment...")
						bridge_start = previous_end_in_genome + 1
						bridge_end = current_start_in_genome - 1
						dna_bridge = extract_fragments_from_scaffold(scaffold, [bridge_start, bridge_end])
						
						# Translate it in the possible reading frames
						modulo_nt = len(dna_bridge) % 3
						if modulo_nt == 0:
							# aa_bridges_list = [dna_bridge.translate()]
							trimmed_dna_bridges = [dna_bridge]
						elif modulo_nt == 1:
							# aa_bridges_list = [dna_bridge[1:].translate, dna_bridge[:-1].translate]
							trimmed_dna_bridges = [dna_bridge[1:], dna_bridge[:-1]]
						elif modulo_nt == 2:
							# aa_bridges_list = [dna_bridge[2:].translate, dna_bridge[:-2].translate]
							trimmed_dna_bridges = [dna_bridge[2:], dna_bridge[:-2]] #for this case it is also possible leaving 1nt at each side, but this leaves 2 gaps instead of one so we exlclude it, it is less parimonious


						# Align each translation and check which one has highest PID
						best_aa_bridge = None
						for i in range(0, len(trimmed_dna_bridges)):
							trimmed_dna_bridge = trimmed_dna_bridges[i]
							aa_bridge = trimmed_dna_bridge.translate()
							print("TRIMMED FRAGMENT (FULL):")
							print(candidate_peptide[fragment_index-1][3])
							print("BRIDGE:")
							print(aa_bridge)

							print("HOMOLOG GAP SEQ:")
							print(homolog_gap_seq.seq)
							
							alignment = align_bridges(homolog_gap_seq, aa_bridge)

							print(alignment[0])
							print(alignment[1])
							
							subs_num = number_of_substitutions(alignment[0], alignment[1])

							print("NUMBER OF SUBSTITUTIONS:")
							print(subs_num)

							aa_bridge_pid = subs_num/homolog_gap

							if best_aa_bridge is None or aa_bridge_pid > best_aa_bridge_pid:
								best_aa_bridge = aa_bridge
								best_aa_bridge_pid = aa_bridge_pid
								best_trimmed_dna_bridge = trimmed_dna_bridge
								best_bridge_index = i
						
						# Correct bridge start and end based on modulo + best bridge index
						if modulo_nt == 1:
							if best_bridge_index == 0:
								corrected_bridge_start = bridge_start + 1
								corrected_bridge_end = bridge_end
							else:
								corrected_bridge_start = bridge_start
								corrected_bridge_end = bridge_end - 1
						elif modulo_nt == 2:
							if best_bridge_index == 0:
								corrected_bridge_start = bridge_start + 2
								corrected_bridge_end = bridge_end
							else:
								corrected_bridge_start = bridge_start
								corrected_bridge_end = bridge_end - 2
						else:
							corrected_bridge_start, corrected_bridge_end = bridge_start, bridge_end

						# Build the entry for the fragments list for the bridge (could be done by extending one of the others also, this seemed easier though)
						if forward:
							best_aa_bridge_entry = [corrected_bridge_start, corrected_bridge_end, best_trimmed_dna_bridge,
							                         best_aa_bridge, previous_end_in_homolog + 1, current_start_in_homolog - 1]
						else:
							best_aa_bridge_entry = [corrected_bridge_end, corrected_bridge_start, best_trimmed_dna_bridge,
							                         best_aa_bridge, previous_end_in_homolog + 1, current_start_in_homolog - 1]

						# Insert new fragment into peptide
						reconstructed_peptides[peptide_index].insert(fragment_index, best_aa_bridge_entry)


						# upstream_superbridge_start = previous_end_in_genome - gaps_difference + 1 # NEED TO ADD BRIDGE START AND END (NOT SUPER) AND GET THE NONSUPER SEQ TO ADD THE X's AND COMPARE
						# upstream_superbridge_end = current_start_in_genome - 1
						# if forward:
						# 	dna_upstream_superbridge = extract_fragments_from_scaffold(scaffold, [upstream_superbridge_start, upstream_superbridge_end])
						# else:
						# 	dna_upstream_superbridge = extract_fragments_from_scaffold(scaffold, [upstream_superbridge_end, upstream_superbridge_start])
						# aa_upstream_superbridge = dna_upstream_superbridge.translate()
						
						# print("UPSTREAM SUPERBRIDGE:")
						# print(aa_upstream_superbridge)
						
						# assert len(aa_upstream_superbridge) == homolog_gap
						# unbridged_sequence = current_fragment[3][:-gaps_difference] + homolog_gap*'X' + candidate_peptide[index-1][3]
						# bridged_sequence = current_fragment[3][:-gaps_difference] + aa_upstream_superbridge + candidate_peptide[index-1][3]
						# assert len(bridged_sequence) == len(unbridged_sequence)

						# subs_num = number_of_substitutions(homolog_gap_seq.seq, aa_upstream_superbridge)
						# upstream_bridged_sequence_pid = subs_num/homolog_gap

						# print("Cutting from next fragment...")
						# downstream_superbridge_start = previous_end_in_genome + 1
						# downstream_superbridge_end = current_start_in_genome - 1 + gaps_difference
						# if forward:
							# dna_downstream_superbridge = extract_fragments_from_scaffold(scaffold, [downstream_superbridge_start, downstream_superbridge_end])
						# else:
							# dna_downstream_superbridge = extract_fragments_from_scaffold(scaffold, [downstream_superbridge_end, downstream_superbridge_start])
						# aa_downstream_superbridge = dna_downstream_superbridge.translate()
						# assert len(aa_downstream_superbridge) == homolog_gap
						# unbridged_sequence = current_fragment[3] + homolog_gap*'X' + candidate_peptide[index-1][3][gaps_difference:]
						# bridged_sequence = current_fragment[3] + aa_downstream_superbridge + candidate_peptide[index-1][3][gaps_difference:]
						# assert len(unbridged_sequence) == len(bridged_sequence)

						# subs_num = number_of_substitutions(homolog_gap_seq, bridged_sequence)
						# downstream_bridged_sequence_pid = subs_num/homolog_gap

						# if upstream_bridged_sequence_pid > downstream_bridged_sequence_pid:
						# 	if upstream_bridged_sequence_pid > unbridged_sequence_pid:
						# 		print("Fragment cut from previous fragment improves the reconstruction. Keeping that fragment.")
						# 		if forward:
						# 			# The fragment previous to the bridge needs to be trimmed, since the 3' end now goes in the superbridge
						# 			reconstructed_peptides[peptide_index][fragment_index-1][1] = upstream_superbridge_start - 1
						# 			reconstructed_peptides[peptide_index][fragment_index-1][2] = reconstructed_peptides[peptide_index][fragment_index-1][2][:-gaps_difference]
						# 			reconstructed_peptides[peptide_index][fragment_index-1][3] = reconstructed_peptides[peptide_index][fragment_index-1][3][:-gaps_difference/3]
						# 			reconstructed_peptides[peptide_index][fragment_index-1][5] = reconstructed_peptides[peptide_index][fragment_index-1][5] - gaps_difference/3

						# 			# Extending the fragment that goes after the gap with the bridge
						# 			reconstructed_peptides[peptide_index][fragment_index][0] = upstream_superbridge_start
						# 			reconstructed_peptides[peptide_index][fragment_index][2] = dna_upstream_superbridge + reconstructed_peptides[peptide_index][fragment_index][2]
						# 			reconstructed_peptides[peptide_index][fragment_index][3] = aa_upstream_superbridge + reconstructed_peptides[peptide_index][fragment_index][3]
						# 			reconstructed_peptides[peptide_index][fragment_index][4] = reconstructed_peptides[peptide_index][fragment_index][4] - homolog_gap

						# 		else:
						# 			reconstructed_peptides[peptide_index][fragment_index-1][0] = upstream_superbridge_start - 1
						# 			reconstructed_peptides[peptide_index][fragment_index-1][2] = reconstructed_peptides[peptide_index][fragment_index-1][2][:-gaps_difference]
						# 			reconstructed_peptides[peptide_index][fragment_index-1][3] = reconstructed_peptides[peptide_index][fragment_index-1][3][gaps_difference/3 - 1:]
						# 			reconstructed_peptides[peptide_index][fragment_index-1][4] = reconstructed_peptides[peptide_index][fragment_index-1][5] + gaps_difference/3

						# 			# Extending the fragment that goes after the gap with the bridge, reverse strand
						# 			reconstructed_peptides[peptide_index][fragment_index][1] = upstream_superbridge_start
						# 			reconstructed_peptides[peptide_index][fragment_index][2] = dna_upstream_superbridge + reconstructed_peptides[peptide_index][fragment_index][2]
						# 			reconstructed_peptides[peptide_index][fragment_index][3] = reconstructed_peptides[peptide_index][fragment_index][3] + aa_upstream_superbridge
						# 			reconstructed_peptides[peptide_index][fragment_index][4] = reconstructed_peptides[peptide_index][fragment_index][4] + homolog_gap

						# 	else:
						# 		print("Not cutting is the best situation.")
						# elif downstream_bridged_sequence_pid > unbridged_sequence_pid:
						# 	print("Fragment cut from next fragment improves the reconstruction. Keeping that fragment.")
						# 	if forward:
						# 		reconstructed_peptides[peptide_index][fragment_index]
						# 	else:
						# 		pass

						# else:
						# 	print("Not cutting is the best situation.")


				else:
					print('No gap in either genome or homolog so no gap can/needs to be filled.')
			else:
				print('Fragment is the first; checking orientation of the strand.')
				if current_fragment[0] < current_fragment[1]:
					print('Fragment is on the forward strand.')
					forward = True
				else:
					print('Fragment is on the reverse strand.')
					forward = False

			# afterwards get the end of the fragment in both the genome and the homolog to account it for the next one
			previous_end_in_homolog = current_fragment[5]
			if forward:
				previous_end_in_genome = current_fragment[1]
			else:
				previous_end_in_genome = current_fragment[0]

	print("FINISHED GAP BRIDGING. OBTAINED PEPTIDES:")
	print(reconstructed_peptides)
	return reconstructed_peptides

def peptides2fasta(reconstructed_peptides):
	"""
	Reconstruct the paleo-peptides from the found fragments, connecting the gaps with X.

	Arguments:
		reconstructed_peptides: list containing the peptides reconstructed from a given seed.

	Returns:
		A list containing the full aminoacid sequences for each of the peptides.
	"""
	print("PUTTING X AT GAPS")
	print(reconstructed_peptides[0])
	peptide_sequences_list = []
	for peptide in reconstructed_peptides:
		full_seq = ''
		previous_end_in_homolog = None
		for fragment in peptide:
			# print(fragment)
			aa_seq = fragment[3]
			start_in_homolog = fragment[4]
			if previous_end_in_homolog is not None:
				gap_size = start_in_homolog - previous_end_in_homolog - 1
			else:
				gap_size = 0
			full_seq = full_seq + gap_size*'X' + aa_seq
			previous_end_in_homolog = fragment[5]
		peptide_sequences_list.append(full_seq)
	return peptide_sequences_list

def blastp(query, target, wordsize, matrix, max_evalue, threads, outprefix, block_size, diamond = False):
	"""
	Run blastp of the reconstructed peptide sequences against a protein database.

	Arguments:
		query: path to FASTA file containing the reconstructed peptide sequences to validate
		target: path to protein database to use for the validation (should be NCBI nr or similar)
		wordsize: word size to use for the blastp run.
		matrix: alignment scoring matrix to use for the blastp run.
		max_evalue: maximum evalue that a hit can have to be shown in the blastp output.
		threads: number of threads to use for the search.
		outprefix: prefix to use for the output file.

	Returns:
		A Pandas dataframe with the results of blastp in tabular format.
	"""
	if diamond:
		blast_file = outprefix + ".diamond_blastp." + matrix + ".evalue" + max_evalue + ".out"
		blast_file = blast_file.replace(' ', '_')
		# blastp_command = "diamond blastp --more-sensitive --max-target-seqs 500 --evalue " + max_evalue + " --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames -b 30 -c 1 --threads " + threads + " -d " + target + " -q " + query + " -o " + blast_file
		blastp_command = "diamond blastp --more-sensitive --max-target-seqs 500 --max-hsps 0 --evalue " + max_evalue + " --outfmt 6 qseqid qlen sallseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames -b " + block_size + " -c 1 --threads " + threads + " -d " + target + " -q " + query + " -o " + blast_file
		# WARNING: CHANGE THE BLAST COMMAND LATER! THE GOOD ONE IS THE COMMETED OUT BUT IN THE LOCAL DB I DO NOT HAVE TAXID INFO FOR THE OUTPUT!
		# NOTE: On Robert's pipeline they use bitscore instead of evalue to filter results: --min-score 40; to keep in mind
		# NOTE 2: We don't have sacc on the output columns since apparantly diamond does not allow for it
	else:
		blast_file = outprefix + ".blastp.wordsize" + wordsize + "." + matrix + ".evalue" + max_evalue + ".out"
		blast_file = blast_file.replace(' ', '_')
		blastp_command="blastp -query " + query + " -db " + target + " -word_size " + wordsize + " -matrix " + matrix + " -evalue " + max_evalue + " -outfmt \"6 qseqid qlen sallseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore stitle staxids sscinames\" -num_threads " + threads + " -out " + blast_file
	if not os.path.exists(blast_file):
		os.system(blastp_command)
	else:
		print("WARNING: BLASTp/Diamond results found at " + blast_file + ". Skipping BLASTp/Diamond step.")
	try:
		blastp_results = pd.read_csv(blast_file, sep='\t', header = None, dtype = {14: 'str'}) # specify str for staxids since there can be multiple ones separated by semicolons
	except pd.errors.EmptyDataError:
		sys.exit("No hits obtained from blastp, exiting program.")
	blastp_results.rename(columns={0: 'qseqid', 1:'qlen', 2: 'sallseqid', 3: 'slen', 4: 'pident', 5: 'length', 6: 'mismatch', 7: 'gapopen', 8: 'qstart', 9: 'qend', 10: 'sstart', 11: 'send', 12: 'evalue', 13: 'bitscore', 14: 'stitle', 15: 'staxids', 16: 'sscinames'}, inplace = True)
	return(blastp_results)

def is_child(query_taxid, parent_taxid, taxdb):
	"""
	Check if a given query taxid is indeed a child of another parent taxid.

	Arguments:
		query_taxid: taxid whose parents nodes we wish to know.
		parent_taxid: taxid for which to check if it is a parent or not.
		taxdb: path to a TaxDb object from taxopy, containing the NCBI Taxonomy database.

	Returns:
		A Boolean value indicating whether or not query_taxid is a child node of parent_taxid.
	"""
	query_taxon = taxopy.Taxon(query_taxid, taxdb)
	parent_taxon = taxopy.Taxon(parent_taxid, taxdb)
	lca_taxon = taxopy.find_lca([query_taxon, parent_taxon], taxdb)
	
	if lca_taxon.taxid == parent_taxon.taxid:
		return True
	else:
		return False

def find_functional(homologs_length_dict, blastp_filtered_summary):
	"""
	Compare the candidate peptides to their homologs to establish if they are likely to be functional or pseudogenised.

	Arguments:
		homologs_length_dict: Dictionary containing the length in AA residues (value) of each protein of the closest homolog's proteome (keys).
		blastp_filtered_summary: Pandas dataframe containing the reconstructed peptides that pass the blastp filtering step (output of filter_blastp_output)
	"""
	# p1 = re.compile('\.pseudopeptide_candidate_[0-9]+')
	# p2 = re.compile('scaffold_[0-9]+\.')
	homologs_length_list = []
	length_ratios_list = []
	for index, row in blastp_filtered_summary.iterrows():
		current_homolog = row['query'].split('___')[1]
		# current_homolog = p1.sub('', p2.sub('', row['query']))
		current_homolog_len = int(homologs_length_dict[current_homolog])
		current_peptide_len = int(row['length (aminoacid)'])
		homologs_length_list.append(current_homolog_len)
		length_ratios_list.append(current_peptide_len/current_homolog_len)
	blastp_filtered_summary['homolog length (aminoacid)'] = homologs_length_list
	blastp_filtered_summary['lengths ratio'] = length_ratios_list
	return blastp_filtered_summary

def filter_blastp_output(blastp_df, parent_taxid, homologs_length_dict, taxdb_nodes = None, taxdb_names = None, taxdb_merged = None, excluded_taxids_list = []):
	"""
	Filter the output of blastp based on the taxonomic assignment of the hits.

	Arguments:
		blastp_df: Pandas dataframe containing the results of BLASTP.
		parent_taxid: taxid of the taxon to which the hits should belong.
		taxdb_nodes: path to the nodes.dmp file of the NCBI Taxonomy database.
		taxdb_names: path to the names.dmp file of the NCBI Taxonomy database
		taxdb_merged: path to the merged.dmp file of the NCBI Taxonomy database.

	Returns: a Pandas dataframe with the results of BLASTP for those queries with at least one hit belonging to the parent_taxid clade, and a second dataframe with a summary of the BLASTP results for those queries.
	"""
	if taxdb_nodes is not None and taxdb_names is not None and taxdb_merged is not None:
		taxdb = taxopy.TaxDb(nodes_dmp = taxdb_nodes, names_dmp = taxdb_names, merged_dmp = taxdb_merged)
	else:
		try:
			taxdb = taxopy.TaxDb
		except taxopy.exceptions.DownloadError:
			sys.exit("ERROR: Failed to download the NCBI Taxonomy database. We recommend downloading it manually and specifying the path to the files.")
	queries = list(set(blastp_df['qseqid']))
	queries.sort()
	peptides_to_keep = []
	alien_indexes = {}
	blastp_summary = pd.DataFrame(data = None, index = [*range(len(queries))], columns = ['query', 'length (aminoacid)', 'belonging_hits_count', 'nonbelonging_hits_count', 
		'belonging_min_eval', 'nonbelonging_min_eval', 'alien_index', 'position_in_scaffold'])
	current_index = 0
	for query in queries:
		print('QUERY:')
		print(query)
		blastp_summary['query'][current_index] = query
		protein_homolog_name = query.split("___")[1]
		print("HOMOLOG NAME:")
		print(protein_homolog_name)
		protein_homolog_seqid = query.split('___')[1].split(".")[1].split("_")[1]
		print("PROTEIN HOMOLOG SEQ ID:")
		print(protein_homolog_seqid)
		scaffold = query.split("___")[0]
		gff_name=".".join(["pseudogene_finder", protein_homolog_name, "reconstructed_peptides.gff"])
		coordinates=""
		with open("extended_peptides_all_gff/"+gff_name, 'r') as fh:
			for gff_entry in fh:
				gff_entry = gff_entry.split("\t")
				if gff_entry[0] == scaffold:
					peptide_number = gff_entry[-1].split(";")[0].strip("ID=pseudogene_")
					print("___".join([scaffold, protein_homolog_name, "pseudopeptide_candidate_" + peptide_number]))
					if query == "___".join([scaffold, protein_homolog_name, "pseudopeptide_candidate_" + peptide_number]):
						print(coordinates)
						if len(coordinates) == 0:
							orientation = gff_entry[6]
							if gff_entry[6] == '+':
								coordinates = "..".join(gff_entry[3:5])
							else:
								coordinates = "c("+"..".join(gff_entry[3:5])
						else:
							assert orientation == gff_entry[6]
							coordinates = coordinates+","+"..".join(gff_entry[3:5])
			if orientation == '-':
				coordinates = coordinates + ")"

		correct_taxa = False
		homolog_in_hits = False
		df_subset = blastp_df.loc[blastp_df['qseqid'] == query]
		belonging_query_min_eval = -1
		nonbelonging_query_min_eval = -1
		belonging_hits_count = 0
		nonbelonging_hits_count = 0
		for index, row in df_subset.iterrows():
			query_taxid = str(row['staxids'])
			query_seqids = str(row['sallseqid'])
			# check if the target protein homolog is among the hits of that particular candidate peptide
			if protein_homolog_seqid in query_seqids:
				homolog_in_hits = True
			if query_taxid not in excluded_taxids_list: # make sure this taxid is not of the ones we want to exclude
				query_hit_evalue = float(row['evalue'])
				if ';' in query_taxid:
					taxids_list = [int(x) for x in query_taxid.split(';') if x not in excluded_taxids_list] # make sure this taxid is not of the ones we want to exclude
					if len(taxids_list) > 1: # make sure we still have some taxids to look for
						print(taxids_list)
						taxa_list = [taxopy.Taxon(x, taxdb) for x in taxids_list]
						query_taxid = taxopy.find_lca(taxa_list, taxdb).taxid
					elif len(taxids_list) == 1:
						query_taxid = taxids_list[0]
					else:
						continue
				if is_child(query_taxid = int(query_taxid), parent_taxid = int(parent_taxid), taxdb = taxdb):
					correct_taxa = True
					belonging_hits_count += 1
					if belonging_query_min_eval == -1 or belonging_query_min_eval > query_hit_evalue:
						belonging_query_min_eval = query_hit_evalue
				else:
					nonbelonging_hits_count += 1
					if nonbelonging_query_min_eval == -1 or nonbelonging_query_min_eval > query_hit_evalue:
						nonbelonging_query_min_eval = query_hit_evalue
		blastp_summary['belonging_hits_count'][current_index] = belonging_hits_count
		blastp_summary['nonbelonging_hits_count'][current_index] = nonbelonging_hits_count
		blastp_summary['belonging_min_eval'][current_index] = belonging_query_min_eval
		blastp_summary['nonbelonging_min_eval'][current_index] = nonbelonging_query_min_eval
		blastp_summary['length (aminoacid)'][current_index] = df_subset['qlen'][df_subset.index[0]] # taking into accound the df_subset has hits all from same query so qlen will be the same for all rows
		blastp_summary['position_in_scaffold'][current_index] = coordinates
		if belonging_query_min_eval == -1:
			alien_indexes[query] = -np.inf # we use -Inf and Inf when Alien Index cannot be computed because either there are no belonging queries or no non-belonging queries
		elif nonbelonging_query_min_eval == -1:
			alien_indexes[query] = np.inf
		elif nonbelonging_query_min_eval == 0 and belonging_query_min_eval == 0:
			alien_indexes[query] = 0
		elif nonbelonging_query_min_eval == 0:
			alien_indexes[query] = -np.inf
		elif belonging_query_min_eval == 0:
			alien_indexes[query] = np.inf
		else:
			alien_indexes[query] = np.log10(nonbelonging_query_min_eval) - np.log10(belonging_query_min_eval)
		blastp_summary['alien_index'][current_index] = alien_indexes[query]
		if correct_taxa and homolog_in_hits and alien_indexes[query] > 0:
			peptides_to_keep.append(query)
		blastp_subset_df = blastp_df.loc[blastp_df['qseqid'].isin(peptides_to_keep)][blastp_df.columns]
		blastp_subset_df = blastp_subset_df.loc[~blastp_subset_df['staxids'].isin(excluded_taxids_list)][blastp_subset_df.columns]
		current_index += 1
	blastp_summary = blastp_summary.loc[blastp_summary['query'].isin(peptides_to_keep)][blastp_summary.columns]
	blastp_summary = find_functional(homologs_length_dict = homologs_length_dict, blastp_filtered_summary = blastp_summary)
	blastp_summary.to_csv('reconstructed_peptides_all.blastp.summary.tsv', sep = '\t', index = False)
	blastp_filtered_summary = blastp_summary[blastp_summary.columns][blastp_summary['alien_index'] > 0]
	# blastp_filtered_summary.to_csv('reconstructed_peptides_all.blastp.filtered.summary.tsv', sep = '\t', index = False)
	return(blastp_subset_df, blastp_filtered_summary)

def subset_fasta(blastp_filtered_summary):
	"""
	Write new FASTA files containing only those peptides that pass the filtering thresholds.

	Arguments:
		blastp_filtered_summary: Pandas dataframe containing the filtered blastp results.

	Returns:
		None. Writes the FASTA files to disk.
	"""
	os.system("mkdir -p filtered_peptides_fasta")
	seqids2keep = list(blastp_filtered_summary['qseqid'])
	files_list = glob.glob('extended_peptides_all_fasta/*.fasta')
	for file in files_list:
		file_string = ''
		with open(file, 'r') as fh:
			for sequence in FastaIO.FastaIterator(fh):
				if sequence.id in seqids2keep:
					file_string = file_string + ">%s\n%s\n" % (sequence.id, sequence.seq)
		if len(file_string) > 0:
			with open(file.replace('.fasta', '.filtered.fasta').replace('extended_peptides_all_fasta', 'filtered_peptides_fasta'), 'w') as wh:
				wh.write(file_string)

def subset_gff(blastp_filtered_summary):
	"""
	Write new GFF files containing only those peptides that pass the filtering thresholds.

	Arguments:
		blastp_filtered_summary: Pandas dataframe containing the fitered blastp results.

	Returns:
		None. Writes the GFF files to disk.
	"""
	os.system("mkdir -p filtered_peptides_gff")
	seqids2keep = list(blastp_filtered_summary['qseqid'])
	files_list = glob.glob('extended_peptides_all_gff/*.gff')
	for file in files_list:
		file_string = ''
		with open(file, 'r') as fh:
			for line in fh:
				# print(line)
				if line.startswith('##'):
					file_string = file_string + line
				else:
					gff_entry = line.strip().split('\t')
					peptide_number = gff_entry[8].split(';')[0].replace('ID=pseudogene_', '')
					peptide_name = gff_entry[0] + file.replace('extended_peptides_all_gff/', '').replace('pseudogene_finder', '').replace('reconstructed_peptides.gff', '') + 'pseudopeptide_candidate_' + peptide_number
					# print(peptide_name)
					if peptide_name in seqids2keep:
						file_string = file_string + line
		if file_string != '##gff-version 3\n':
			with open(file.replace('.gff', '.filtered.gff').replace('extended_peptides_all_gff', 'filtered_peptides_gff'), 'w') as wh:
				wh.write(file_string)

def check_stop_codons(blastp_results, outprefix):
	"""
	Given the summary of candidate peptides that pass the filtering thresholds, check if there are stop codons in their sequence to assess how likely they are to be pseudogenes.

	Arguments:
		blastp_results: Pandas dataframe containing the summary of the filtered blastp results.
		outprefix: Prefix used to name the output files.

	Returns:
		The same Pandas dataframe with an extra column indicating if the sequences contain or not stop codons.
	"""
	# p1 = re.compile('\.pseudopeptide_candidate_[0-9]+')
	# p2 = re.compile('scaffold_[0-9]+\.')
	has_stop_codons = []
	for index, row in blastp_results.iterrows():
		peptide_id = row['query']
		print("PEPTIDE ID:")
		print(peptide_id)
		homolog = peptide_id.split('___')[1]
		# homolog = p1.sub('', p2.sub('', peptide_id))
		fasta_file = 'filtered_peptides_fasta/%s.%s.reconstructed_peptides.filtered.fasta' % (outprefix, homolog)
		count = 0
		with open(fasta_file) as fh:
			for peptide in FastaIO.FastaIterator(fh):
				if peptide.id == peptide_id:
					count += 1
					if '*' in peptide.seq:
						has_stop_codons.append('Y')
					else:
						has_stop_codons.append('N')
		assert count == 1
	blastp_results['Stop codons presence'] = has_stop_codons
	return(blastp_results)



if __name__ == '__main__':
	__version__ = "0.5.2"
	
	#### PARSE ARGS ####
	args = docopt(__doc__)
	try:
		args['--excluded_taxids'] = [x for x in args['--excluded_taxids'].split(',')]
	except ValueError:
		sys.exit('ERROR: Please make sure --parent_taxid is an integer.')
	except AttributeError:
		pass
	try:
		args['--parent_taxid'] = int(args['--parent_taxid'])
	except TypeError:
		pass
	print(args)


	#### MAIN ####
	#### SUBCOMMAND 1: tblastn ####
	if args['runall'] or args['tblastn']:
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

	if args['tblastn']:
		if args ['--verbose']:
			print("tBLASTn completed.\nExecution finished.")
		exit()

	if args['extend']:
		tblastn_output = pd.read_csv(args['--tblastn_output'], sep='\t', header = None)
		tblastn_output.rename(columns={0: 'qseqid', 1: 'sseqid', 2: 'pident', 3: 'length', 4: 'mismatch', 5: 'gapopen', 6: 'qstart', 7: 'qend', 8: 'sstart', 9: 'send', 10: 'evalue', 11: 'bitscore'}, inplace = True)

	if args['extend'] or args['runall']:
		if args['--verbose']:
			print("Extending the peptide seeds found...")
		# The rest of steps are done protein homolog by protein homolog
		homologs_length_dict = {}
		with open(args['--proteins']) as proteins_fh:
			for protein in FastaIO.FastaIterator(proteins_fh):
				protein_id = protein.id
				print('STARTING NEW PROTEIN: %s' % str(protein_id))
				if '[' in protein_id or ']' in protein_id or '=' in protein_id or '(' in protein_id or ')' in protein_id:
					print("WARNING: protein name contains special characters. They have been replaced. Please make sure this is not a problem and if so change your protein IDs manually")
					protein_id = protein_id.replace('[', '')
					protein_id = protein_id.replace(']', '')
					protein_id = protein_id.replace('(', '')
					protein_id = protein_id.replace(')', '')
					protein_id = protein_id.replace('=', '_')
					protein_id = protein_id.replace('/', '_')

				homologs_length_dict[protein_id] = len(protein)
				if args['--verbose']:
					print("Processing results for protein %s" % protein_id)

				# Get which positions in the query genome that match the protein homolog
				primary_seeds = find_primary_seed_coordinates(tblastn_output, protein.id)
				
				# print(primary_seeds)

				# Extract the nucleotide sequence corresponding to those positions
				extract_fragments_from_fasta(file=args['--genome'], fragments_list=primary_seeds)

				# print(primary_seeds)

				# Translate the nucleotide sequence to AA
				translate_dna(fragments_list=primary_seeds, frame = 0) # since this comes from tblastn the frist frame is already the good one

				# print(primary_seeds)

				# Align the peptide seed to the protein homolog
				# primary_seed_alignment_list = align_peptides(protein = protein, fragments_list=primary_seeds)
				# AlignIO.write(primary_seed_alignment_list, 'alignments.' + protein.id + '.all.fa', 'fasta')

				# Conduct seed extension
				os.system('mkdir -p alignments && mkdir -p extended_peptides_all_fasta && mkdir -p extended_peptides_all_gff')
				with open('extended_peptides_all_gff/' + args['--outprefix'] + '.' + protein_id + '.reconstructed_peptides.gff', 'w') as fh:
					fh.write('##gff-version 3\n')
					with open('extended_peptides_all_fasta/' + args['--outprefix'] + '.' + protein_id + '.reconstructed_peptides.fasta', 'w') as fh_seqs:
						for scaffold, scaffold_candidate_peptides in primary_seeds.items():
							id_num = 0 # id to use later to track the number of peptides per homologous protein and scaffold in the GFF files
							# fh.write(scaffold + '\n')
							print('STARTING NEW SCAFFOLD: ' + scaffold)
							print('RESETTING PEPTIDE COUNT')
							scaffold_seq = get_scaffold_from_fasta(args['--genome'], scaffold)
							# print(len(scaffold_seq))
							for candidate_peptide in scaffold_candidate_peptides: # here candidate peptides means seeds
								print("STARTING NEW SEED:")
								print(candidate_peptide)
								reconstructed_peptides_downstream = []
								reconstructed_peptides_upstream = []
								reconstructed_peptides_downstream.append(candidate_peptide)
								reconstructed_peptides_upstream.append([])
								peptide_index = scaffold_candidate_peptides.index(candidate_peptide)
								# print('CANDIDATE PEPTIDE SEED: ')
								# print(candidate_peptide)
								# print('EXTENDING DOWNSTREAM OF THE GENOME:')

								for peptide_tuple in extend_candidate_peptide(candidate_peptide[0], scaffold_seq, protein, direction = 'downstream', order = 0):
									# print('Yielded peptide tuple is: ')
									# print(peptide_tuple)
									# check the order and based on the length of the already built peptide duplicate or not
									new_peptide = peptide_tuple[0]
									order = peptide_tuple[1]
									# print('downstream peptides:')
									# print(reconstructed_peptides_downstream)
									# print('Length of last downstream peptide: ' + str(len(reconstructed_peptides_downstream[-1])))
									if len(reconstructed_peptides_downstream[-1]) == order:
										# print('Adding new fragment at the end of last peptide.')
										reconstructed_peptides_downstream[-1].append(new_peptide)
										# print('downstream peptides now like this: ')
										# print(reconstructed_peptides_downstream)
										# print('upstream peptides now like this: ')
										# print(reconstructed_peptides_upstream)

									elif len(reconstructed_peptides_downstream[-1]) > order:
										# print('Adding new peptide with the same fragments up to order of new one.')
										reconstructed_peptides_downstream.append(reconstructed_peptides_downstream[-1][0:order] + [new_peptide])
										# print('downstream peptides now like this: ')
										# print(reconstructed_peptides_downstream)
										# print('upstream peptides now like this: ')
										# print(reconstructed_peptides_upstream)
								os.system('mv alignments/pairwise_seqs.tmp.fa alignments/' + protein_id + '_' + scaffold + '_candidate_peptide_' + str(peptide_index) + '_downstream.aln.fa')
								# print('EXTENDING UPSTREAM OF THE GENOME:')
								
								for peptide_tuple in extend_candidate_peptide(candidate_peptide[0], scaffold_seq, protein, direction = 'upstream', order = 0):
									# print('Yielded peptide tuple is:')
									# print(peptide_tuple)
									# check the order and based on the length of the already built peptide duplicate or not
									new_peptide = peptide_tuple[0]
									order = peptide_tuple[1]
									if len(reconstructed_peptides_upstream[-1]) == order - 1:
										reconstructed_peptides_upstream[-1].append(new_peptide)
									elif len(reconstructed_peptides_upstream[-1]) > order - 1:
										reconstructed_peptides_upstream.append(reconstructed_peptides_upstream[-1][0:order - 1] + [new_peptide])
								
								os.system('mv alignments/pairwise_seqs.tmp.fa alignments/' + protein_id + '_' + scaffold + '_candidate_peptide_' + str(peptide_index) + '_upstream.aln.fa')

								# Remove initial seed from upstream part so as not to have it duplicated and turn around the peptide as it was built upside down
								for index in range(0, len(reconstructed_peptides_upstream)):
									reconstructed_peptides_upstream[index].reverse()

								# Combine the upstream and downstream parts of the peptides together
								reconstructed_peptides_complete = []
								print('RECONSTRUCTED PEPTIDES (FRAGMENTED): ')
								print(reconstructed_peptides_upstream)
								print(reconstructed_peptides_downstream)
								for upstream_peptide in reconstructed_peptides_upstream:
									for downstream_peptide in reconstructed_peptides_downstream:
										if upstream_peptide is not None:
											# print('UPSTREAM PEPTIDE IS:')
											# print(upstream_peptide)
											reconstructed_peptides_complete.append(upstream_peptide + downstream_peptide)
										else:
											reconstructed_peptides_complete.append(downstream_peptide)

								# For those that come from the - strand, reverse the order of the fragments again to get the real direction of the peptide
								for index in range(0, len(reconstructed_peptides_complete)):
									this_peptide = reconstructed_peptides_complete[index]
									if this_peptide[0][0] > this_peptide[0][1]:
										reconstructed_peptides_complete[index].reverse()
								# print(reconstructed_peptides_complete)
								
								# Try to fill gaps in those peptides that have some NTs in the genome and also some AAs in the homolog between them
								if not args['--no_gap_bridging']:
									reconstructed_peptides_complete = gap_bridging(reconstructed_peptides_complete, scaffold_seq, homolog = protein)

								# Write reconstructed peptides to a GFF file
								id_num_gff = id_num # rest this to the same value as id_num
								for peptide in reconstructed_peptides_complete:
									fragment_num = 0
									for fragment in peptide:
										fragment_num += 1
										if fragment[0] < fragment[1]:
											genomic_start = fragment[0]
											genomic_end = fragment[1]
											strand = '+'
										else:
											genomic_start = fragment[1]
											genomic_end = fragment[0]
											strand = '-'
										attributes = 'ID=pseudogene_%s;fragment=fragment_%s' % (str(id_num_gff), str(fragment_num))
										gff_string = '\t'.join([scaffold, '.', 'pseudogene', str(genomic_start), str(genomic_end), '.', 
											strand, '.', attributes])
										fh.write(str(gff_string) + '\n')
									id_num_gff += 1 # since there can be more than one peptide per seed, we also increase here for each peptide reconstructed from a given seed, but we increase locally because we need to keep the previous count for the FASTA output
								
								# Write reconstructed peptide sequences to FASTA file:
								peptide_sequences = peptides2fasta(reconstructed_peptides_complete)	
								# print('RECONSTRUCTED PEPTIDES (JOINED SEQUENCES):')
								# print(peptide_sequences)
								# print('NUMBER OF RECONSTRUCTED PEPTIDE SEQUENCES:')
								# print(len(peptide_sequences))
								for index in range(0, len(peptide_sequences)):
									peptide = peptide_sequences[index]
									line = '>%s___%s___pseudopeptide_candidate_%s\n%s\n' % (scaffold, protein_id, str(index+id_num), peptide)
									fh_seqs.write(line)
								# now this will start a new seed potentially in the same scaffold and homologous peptide, so we increase by one just in case
								# print('INCREASING PEPTIDE COUNT BY: %s' % str(len(peptide_sequences)))
								id_num += len(peptide_sequences)
		# Remove temporary files used for alignments
		os.system('rm *.fa && rm lalign.aln && rm lalign.log')
	
	if args['extend']:
		if args['--verbose']:
			print('Seed extension phase completed.\nExecution finished.')
		sys.exit()

	if args['blastp'] or args['runall']:
		if args['--verbose']:
			print("Conducting BLASTp to validate the reconstructed peptides...")
		# Run blastp validation of the reconstructed peptide sequences (putting them on a single file for speed)
		os.system('rm -f reconstructed_peptides_all.fasta && for file in extended_peptides_all_fasta/*.reconstructed_peptides.fasta; do cat $file >> reconstructed_peptides_all.fasta; done')
		if args['--diamond']:
			blastp_output = blastp(query = 'reconstructed_peptides_all.fasta', target = args['--blastp_db'], wordsize = args['--blastp_wordsize'], matrix = args['--blastp_matrix'],
				max_evalue = args['--blastp_max_evalue'], threads = args['--blastp_threads'], outprefix = 'reconstructed_peptides_all', block_size = args['--diamond_block_size'], diamond = True)
		else:
			blastp_output = blastp(query = 'reconstructed_peptides_all.fasta', target = args['--blastp_db'], wordsize = args['--blastp_wordsize'], matrix = args['--blastp_matrix'],
				max_evalue = args['--blastp_max_evalue'], threads = args['--blastp_threads'], outprefix = 'reconstructed_peptides_all', block_size = args['--diamond_block_size'])

		# Remove file with all sequences concatenated which is not needed
		os.system('rm -f reconstructed_peptides_all.fasta')

	if args['blastp']:
		if args['--verbose']:
			print("BLASTp completed.\nExecution finished.")
		sys.exit()

	if args['filter_blastp']:
		blastp_output = pd.read_csv(args['--blastp_output'], sep='\t', header = None, dtype = {14: 'str'}) # specify str for staxids because there can be more than one seprated by semicolons
		blastp_output.rename(columns={0: 'qseqid', 1: 'qlen', 2: 'sallseqid', 3: 'slen', 4: 'pident', 5: 'length', 6: 'mismatch', 7: 'gapopen', 8: 'qstart', 9: 'qend', 10: 'sstart', 11: 'send', 12: 'evalue', 13: 'bitscore', 14: 'stitle', 15: 'staxids', 16: 'sscinames'}, inplace = True)
		blast_file = args['--blastp_output'].replace('.out', '.filtered_taxid' + str(args['--parent_taxid']) + '.out')
		homologs_length_dict = {}
		with open(args['--proteins']) as proteins_fh:
			for protein in FastaIO.FastaIterator(proteins_fh):
				protein_id = protein.id
				if '[' in protein_id or ']' in protein_id or '=' in protein_id or '(' in protein_id or ')' in protein_id:
					print("WARNING: protein name contains special characters. They have been replaced. Please make sure this is not a problem and if so change your protein IDs manually")
					protein_id = protein_id.replace('[', '')
					protein_id = protein_id.replace(']', '')
					protein_id = protein_id.replace('(', '')
					protein_id = protein_id.replace(')', '')
					protein_id = protein_id.replace('=', '_')
					protein_id = protein_id.replace('/', '_')
				homologs_length_dict[protein_id] = len(protein)

	elif args['runall'] and args['--diamond']:
		blast_file = "reconstructed_peptides_all" + ".diamond_blastp." + str(args['--blastp_matrix']) + ".evalue" + str(args['--blastp_max_evalue']) + ".filtered_taxid" + str(args['--parent_taxid']) + ".out"
	elif args['runall'] and not args['--diamond']:
		blast_file = "reconstructed_peptides_all" + ".blastp.wordsize" + str(args['--blastp_wordsize']) + "." + str(args['--blastp_matrix']) + ".evalue" + str(args['--blastp_max_evalue']) + ".filtered_taxid" + str(args['--parent_taxid']) + ".out"

	if args['filter_blastp'] or args['runall']:		
		if args['--parent_taxid'] != 1:
			if args['--verbose']:
				print("Filtering BLASTp output...")
			if args['--excluded_taxids'] is not None:
				blastp_results, blastp_filtered_summary = filter_blastp_output(blastp_output, args['--parent_taxid'], homologs_length_dict, os.path.dirname(__file__).strip('.') + '/nodes.dmp', os.path.dirname(__file__).strip('.') + '/names.dmp', os.path.dirname(__file__).strip('.') + '/merged.dmp', excluded_taxids_list = args['--excluded_taxids'])
			else:
				blastp_results, blastp_filtered_summary = filter_blastp_output(blastp_output, args['--parent_taxid'], homologs_length_dict, os.path.dirname(__file__).strip('.') + '/nodes.dmp', os.path.dirname(__file__).strip('.') + '/names.dmp', os.path.dirname(__file__).strip('.') + '/merged.dmp')
			blastp_results.to_csv(blast_file, sep = '\t', index = False)
			subset_fasta(blastp_results)
			subset_gff(blastp_results)
			blastp_filtered_summary = check_stop_codons(blastp_filtered_summary, args['--outprefix'])
			blastp_filtered_summary.to_csv('reconstructed_peptides_all.blastp.filtered.summary.tsv', sep = '\t', index = False)

		if args['--verbose']:
			print("Filtering of BLASTP results completed.\nExecution finished.")
			sys.exit()

#### END ####
