#!/usr/bin/env python3

#This file is part of Mutational Signatures Project.

#Mutational Signatures Project: need info on project

#Copyright (C) 2018 Erik Bergstrom

#

#Mutational Signatures is free software: need to include distribtution

#rights and so forth

#

#Mutational Signatures is distributed in the hope that it will be useful,

#but WITHOUT ANY WARRANTY; without even the implied warranty of

#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the

#GNU General Public License for more details [change to whatever group we should include.
 

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu

import sys
import os
import re

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])


def context_distribution (context_input, output_file, chromosome_path, chromosome_TSB_path, chromosomes):
	'''
	Creates a csv file for the distribution of nucleotides given a specific context. 
	This csv file needs to be created before simulating mutationalsigantures for the 
	given context.

	Requires: 
		Chromosomes saved in individual text files ('X.txt','Y.txt', '1.txt', etc)
		Transcriptional data saved in binary files for each chromosome. These files
			can be created using the script: "save_tsb_192.py"

	Parameters:
			  context_input  -> simulation context of interest (ex: 96, 192, 1536, 3072, DINUC, INDEL)
				output_file  -> file where the distribution for the given nucleotide context is saved (csv file)
			chromosome_path  -> path to the reference chromosomes
		chromosome_TSB_path  -> path to the transcriptional strand bias references
				chromosomes  -> list of chromosomes for the species of interest
	Returns:
		None

	Outputs:
		CSV file with each chromosome represented as a column and reach row 
		represented as a nucleotide. Under each chromosome for a given nucleotide
		is the proportion of the total length of that chromosome associated with that 
		nucleotide. 
	'''


	dinuc_types = ['AA', 'AC', 'AG', 'AT', 'CA', 'CC', 'CG', 'GA', 'GC', 'TA']


	# Set the context parameter based upon the user input
	if context_input == '96' or context_input == '192':
		context = 3
	elif context_input == '1536':
		context = 5
	elif context_input == 'DINUC':
		context = 2
	else:
		print('Not a valid context')
		sys.exit()


	probs = {}
	chromosome_lengths = []

	# Iterate through each chromosome and open the associated file
	for chrom in chromosomes:
		with open (chromosome_path + chrom + ".txt") as f:
			chromosome = f.readline().strip()
			chromosome_lengths.append(len(chromosome))

			# If desired context is for TSB, open the transcriptional file, too
			if context_input == '192':
				with open (chromosome_TSB_path + chrom + "_192.txt", 'rb') as tsb:
					chromosome_tsb = tsb.read()

			# Iterate through the chromosome base by base
			for i in range (0, (len(chromosome)-context), 1):
				nuc = chromosome[i:i+context]
				base = nuc[int(context/2)]

				# Skip the base if unknown 
				if "N" in nuc:
					pass

				else:
					if context_input != "DINUC":
						# Only save the pyrimidine context (canonical)
						if base == 'A' or base == 'G':
							nuc = revcompl(nuc)

						# Adjust the nucleotide representaiton if TSB is desired
						if context_input == '192':
							bias = chromosome_tsb[i+int(context/2)]

							if bias == 0:
								nuc = 'N:' + nuc
							elif bias == 1:
								nuc = 'T:' + nuc
							elif bias == 2:
								nuc = 'U:' + nuc
							else:
								nuc = 'B:' + nuc 

						# Update the dictionary for the current nucleotide
						if nuc not in probs.keys():
							probs[nuc] = {chrom:1}
						else:
							if chrom not in probs[nuc].keys():
								probs[nuc][chrom] = 1
							else:
								probs[nuc][chrom] += 1

					else:
						# Populate the dictionary if desired context is DINUC
						for dinuc in dinuc_types:
							probs[dinuc] = {}

						# Update the dictionary for the current nucleotide
						if nuc not in probs.keys():
							nuc = revcompl(nuc)
							
						if chrom not in probs[nuc]:
							probs[nuc][chrom] = 1
						else:
							probs[nuc][chrom] += 1


		print("chrom ", chrom, "done")

	# Write the resulting dictionary to the csv file
	with open (output_file, 'w') as out:
		print(' ,', end='', file=out)
		for chrom in chromosomes[:-1]:
			print(chrom + ',', end='',file=out)
		print(chromosomes[-1],file=out)
		for nuc in probs.keys():
			print (nuc + ',', end='', file=out)
			for i in range (0, len(chromosomes[:-1]), 1):
				print(str(probs[nuc][chromosomes[i]]/chromosome_lengths[i]) + ',', end='', file=out)
				out.flush()
			print(probs[nuc][chromosomes[-1]]/chromosome_lengths[-1], file=out)


	# Sort the file so that the nucleotides are in alphabetical order
	sort_command_1 = "sort -t ',' -k 1,1 "
	sort_command_2 = " -o "
	os.system (sort_command_1 + output_file + sort_command_2 + output_file)




def main():
	context_input = sys.argv[1]
	genome = sys.argv[2]


	chromosomes = ['X', 'Y', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
			   	   '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']

	if genome.upper() == 'MM10':
		chromosomes = chromosomes[:19]

	script_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', script_dir)
	chromosome_path = ref_dir + "/references/chromosomes/chrom_string/" + genome + "/"
	chromosome_TSB_path = ref_dir + "/references/chromosomes/tsb/" + genome + "/"

	output_file = ref_dir + '/references/chromosomes/context_distributions/context_distribution_' + genome + "_" + context_input + '.csv'
	context_distribution(context_input, output_file, chromosome_path, chromosome_TSB_path, chromosomes)
    
if __name__ == '__main__':
    main()

