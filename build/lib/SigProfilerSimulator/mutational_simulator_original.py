#!/usr/bin/env python3

#Author: Erik Bergstrom

#Contact: ebergstr@eng.ucsd.edu


import time
import sys
import random
import fastrand
import os
import pickle
import subprocess
import argparse
import datetime
import shutil
import bisect
import numpy as np
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as matRef
# from memory_profiler import profile
# from pympler.tracker import SummaryTracker
import pandas as pd
import re
start_run = time.time()


#################################### Functions ###########################################


revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N'}[B] for B in x][::-1])
revbias = lambda x: ''.join([{'0':'0', '3':'3', '1':'2','2':'1','U':'T','T':'U','B':'B','N':'N','Q':'Q'}[B] for B in x][::-1])


def noise (samples, noisePoisson=False, noiseAWGN=0):
	if noisePoisson:
		for mut in samples:
			noise_value = np.random.poisson(samples[mut])
			samples[mut] = noise_value
		return(samples)

	if noiseAWGN != 0:
		for mut in samples:
			lower_bound = samples[mut] - int(noiseAWGN/2)
			upper_bound = samples[mut] + int(noiseAWGN/2)
			noise_value = int(np.random.uniform(lower_bound, upper_bound))
			if noise_value < 0:
				noise_value = 0
			samples[mut] = noise_value
		return(samples)


# @profile
def combine_simulation_files (iterations, output_path, chromosomes, samples=[], bed=False, exome=False, vcf=False):
	'''
	Combines the separate sample_iteration_chromosome simulated files into a 
	single file per sample per iteration.

	Parameters:
		 iterations  -> Number of iterations that were simulated per sample
		output_path  -> The path where the simulations are located
		chromosomes  -> A list of the chromosomes that were simulated
		    samples  -> A list of the sample names that were simulated

	Returns:
		None
	Outputs:
		-> single vcf files per sample per iteration
	'''

	if vcf:
		for sample in samples:
			for i in iterations:
				with open(output_path + sample + "/" + sample + "_" + str(i) + ".vcf", "wb") as f:
					for chrom in chromosomes:
						if True:


						# try:
							with open(output_path + sample + "/" + sample + "_" + str(i) + "_" + chrom + ".vcf", "rb") as fd:
								shutil.copyfileobj(fd, f)
							os.remove(output_path + sample + "/" + sample + "_" + str(i) + "_" + chrom + ".vcf")
						# except:
						# 	pass


	else:		
		extension = ''
		if bed and not exome:
			extension = "_BED"
		for i in iterations:
			with open(output_path + str(i) + ".maf", "wb") as f:
				for chrom in chromosomes:
					try:
						with open(output_path + str(i) + "_" + chrom + extension + ".maf",'rb') as fd:
							shutil.copyfileobj(fd, f)
						os.remove(output_path  + str(i) + "_" + chrom + extension + ".maf")
					except:
						pass

def chrom_proportions (chrom_path, genome, chromosomes):
	'''
	Creates a text file that contains the proportional size of each chromosome in relation to 
	the entire genome. The proportions are saved into a list.

	Parameters:
		 chrom_path  -> path to the chromosome string files
			 genome  -> name of the genome of interest
		chromosomes  -> list of chromosomes for the species of interest

	Returns:
		None

	Outputs:
		-> a text file saved into the chrom_path with the name:
		   genome + _proportions.txt (ex: GRCh37_proportions.txt)
	'''
	chromosome_lengths = []
	chromosomeProbs = []
	total_length = 0
	for chrom in chromosomes:
		with open (chrom_path + chrom + ".txt", "rb") as f:
			chromosome = f.read().strip()
			chromosome_lengths.append(len(chromosome))
			total_length += len(chromosome)

	for lengths in chromosome_lengths:
		chromosomeProbs.append(lengths/total_length)

	with open (chrom_path + genome + "_proportions.txt", 'wb') as out:
		pickle.dump(chromosomeProbs, out)

def chrom_proportions_BED (bed_file, chrom_path, genome, chromosomes):
	'''
	Creates a text file that contains the proportional size of each chromosome in relation to 
	the given BED file ranges. The proportions are saved into a list.

	Parameters:
		   bed_file  -> input file that contains the desired ranges
		  chrom_path -> path to the chromosome string files
			 genome  -> name of the genome of interest
		chromosomes  -> list of chromosomes for the species of interest

	Returns:
		None

	Outputs:
		-> a text file saved into the chrom_path with the name:
		   genome + _proportions.txt (ex: GRCh37_proportions.txt)
	'''
	chromosome_lengths = {}
	chromosomeProbs = []
	total_length = 0
	first_line = True
	length = 0
	total = 0
	chromosomes = []
	with open (bed_file) as f:
		next(f)
		for lines in f:
			line = lines.strip().split()
			chrom = line[0]
			if chrom not in chromosomes:
				chromosomes.append(chrom)
			if len(chrom) > 2:
				chrom = chrom[3:]
			start = int(line[1])
			end = int(line[2])
			if first_line:
				chrom_initial = chrom
				first_line = False


			if chrom == chrom_initial:
				length += end-start
			else:
				total += length
				chromosome_lengths[chrom_initial] = length
				chrom_initial = chrom
				length = end-start
				total_length += length
		total_length += length
		chromosome_lengths[chrom_initial] = length

	for chroms in chromosomes:
		chromosomeProbs.append(chromosome_lengths[chroms]/total_length)

	with open (chrom_path + "BED_" + genome + "_proportions.txt", 'wb') as out:
		pickle.dump(chromosomeProbs, out)

# @profile
def update_chromosome ( chrom, location, bases, context):
	'''
	Updates a given chromosome or sequence based upon a given context.
	
	Parameters:
		   chrom -> sequence
		location -> starting position of desired update in the sequence 
		   bases -> desired bases to update (del, ins, SNP, Dinuc, etc.)
		 context -> simulation context (INDEL, DINUC, SNP)
	
	Returns:
		returns -> updated chromosome
		
	Example:
		
		update_chromosome (1.txt, 10546, 'ACG', 'Ins')
		output -> original chromosome ( ...GAAATCT...) becomes ( ...GAAA[ACG]TCT...)
	'''
	chromosome = chrom

	if context == 'Del':
		for i in range (0, len(bases), 1):
			chromosome[location + i] = int(bases[i])
	elif context == 'Ins':
		for i in range (0, len(bases), 1):
			chromosome[location+1+i] = int(bases[i])
	elif context == 'SNP':
		chromosome[location] = int(bases)

	else:
		for i in range(0, len(bases), 1):
			chromosome[location+i] = bases[i]
	return(chromosome)

# @profile
def random_base (limit_start, limit_end):
	'''
	Returns a random nucleotide.

	Inputs: None

	Returns: A, C, G, or T
	'''

	return (('ATCG')[random.randint(limit_start,limit_end)])

# @profile
def bed_ranges (chrom_sim, bed_file, cushion):
	'''
	Returns a list containing the positions corresponding with the desired
	ranges from the BED file provided by the user.

	Parameters:
		  chrom_sim  -> chromosome of interest
		   bed_file  -> input file that contains the desired ranges

	Returns:
		chrom_range  -> a list with all positions that fall within the 
						desired ranges for a given chromosome.
	'''

	chrom_range = []
	chrom_reached = False
	first_reached = True
	with open(bed_file) as f:
		next(f)
		for lines in f:

			# Saves the relevant data for each line in the file
			line = lines.strip().split()
			chrom = line[0]
			if len(chrom) > 2:
				chrom = chrom[3:]
			start = int(line[1])
			end = int(line[2])

			# Skips the line if it does not match the chromosome
			# of interest
			if chrom != chrom_sim:
				if first_reached:
					try:
						next(f)
					except:
						break
				else:
					break

			# Saves the positions within the given range if the range
			# falls on the chromosome of interest
			else:
				if first_reached:
					chrom_reached = True
					first_reached = False
				for i in range (start-cushion, end+cushion, 1):
					chrom_range.append(i)

	# Sorts and returns the list in case the BED file was not sorted previously
	chrom_range.sort()
	chrom_range = list(set(chrom_range))
	return(chrom_range)


def mutation_preparation_chromosomes (catalogue_files, matrix_path, chromosomes, project, log_file):
	out = open(log_file, 'a')
	sample_names = []
	samples = dict()
	mutation_tracker = {}
	for context in catalogue_files:
		mutation_tracker[context] = {}
		with open (catalogue_files[context]) as f:
			first_line = f.readline().strip().split('\t')
		sample_names += first_line[1:]
		samples[context] = {}
		current_samples = first_line[1:]

		# Save the mutation counts for each sample for each nucleotide context
		with open (catalogue_files[context]) as f:
			next(f)
			for lines in f:
				line = lines.strip().split()
				nuc = line[0]
				if nuc == 'complex' or nuc == 'non-matching':
					continue

				sample_index = 1
				for sample in current_samples:
					mutCount = int(line[sample_index])
					if sample not in samples[context]:
						samples[context][sample] = {nuc:mutCount}
					else:
						samples[context][sample][nuc] = int(mutCount)
					sample_index += 1  
		for chrom in chromosomes:
			with open(catalogue_files[context] + ".chr" + chrom) as f:
				next(f)
				for lines in f:
					line = lines.strip().split()
					nuc = line[0]
					if nuc == 'complex' or nuc == 'non-matching':
						continue

					sample_index = 1
					for sample in current_samples:
						mutCount = int(line[sample_index])
						if sample not in mutation_tracker[context]:
							mutation_tracker[context][sample] = {nuc:{}}
							for chromo in chromosomes:
								mutation_tracker[context][sample][nuc][chromo] = 0
							mutation_tracker[context][sample][nuc][chrom] = int(mutCount)
						else:
							if nuc not in mutation_tracker[context][sample]:
								mutation_tracker[context][sample][nuc] = {}
								for chromo in chromosomes:
									mutation_tracker[context][sample][nuc][chromo] = 0

							mutation_tracker[context][sample][nuc][chrom] = int(mutCount)
						sample_index += 1  


	sample_names = list(set(sample_names))
	print("Files successfully read and mutations collected. Mutation assignment starting now.", file=out)
	print("Files successfully read and mutations collected. Mutation assignment starting now.")
	out.close()
	return (sample_names, samples, mutation_tracker)

def mutation_preparation_region(catalogue_files, matrix_path, project, log_file, region):
	out = open(log_file, 'a')
	sample_names = []
	samples = dict()
	mutation_tracker = {}
	for context in catalogue_files:
		mutation_tracker[context] = {}
		with open (catalogue_files[context]) as f:
			first_line = f.readline().strip().split('\t')
		sample_names += first_line[1:]
		samples[context] = {}
		current_samples = first_line[1:]

		# Save the mutation counts for each sample for each nucleotide context
		with open (catalogue_files[context]) as f:
			next(f)
			for lines in f:
				line = lines.strip().split()
				nuc = line[0]
				if nuc == 'complex' or nuc == 'non-matching':
					continue

				sample_index = 1
				for sample in current_samples:
					mutCount = int(line[sample_index])
					if sample not in samples[context]:
						samples[context][sample] = {nuc:mutCount}
					else:
						samples[context][sample][nuc] = int(mutCount)
					sample_index += 1  
		with open(catalogue_files[context]) as f:
			next(f)
			for lines in f:
				line = lines.strip().split()
				nuc = line[0]
				if nuc == 'complex' or nuc == 'non-matching':
					continue

				sample_index = 1
				for sample in current_samples:
					mutCount = int(line[sample_index])
					if sample not in mutation_tracker[context]:
						mutation_tracker[context][sample] = {nuc:{}}
						# for chromo in chromosomes:
							# mutation_tracker[context][sample][nuc][chromo] = 0
						mutation_tracker[context][sample][nuc][region] = int(mutCount)
					else:
						if nuc not in mutation_tracker[context][sample]:
							mutation_tracker[context][sample][nuc] = {}
							# for chromo in chromosomes:
								# mutation_tracker[context][sample][nuc][chromo] = 0

						mutation_tracker[context][sample][nuc][region] = int(mutCount)
					sample_index += 1  


	sample_names = list(set(sample_names))
	print("Files successfully read and mutations collected. Mutation assignment starting now.", file=out)
	print("Files successfully read and mutations collected. Mutation assignment starting now.")
	out.close()
	return (sample_names, samples, mutation_tracker)


def mutation_preparation (catalogue_files, log_file):
	'''
	Returns a list of all sample names and a dictionary containing the mutation count
		for each mutation type for a given context.
		
	Parameters:
		catalogue_files -> all mutational matrix catalgoues that include the number of mutations
						   for a given mutation type for a given sample

	Returns:
		sample_names  -> list of all sample names
			 samples  -> dictionary containing the mutation count for each mutation type
						 for a given context.

	Example return value:
		sample_names = ['PDXXXX', 'PDYYYY', ...]
		samples = {'96':{'PDXXXX':{'A[A>C]A':35, 'A[A>T]A':12, ...}
						 'PDYYYY':{'A[A>C]A':23, 'A[A>T]A':9,  ...}},
				   'DINUC':{'PDXXXX':{'A[A>C]A':35, 'A[A>T]A':12, ...}
							'PDYYYY':{'A[A>C]A':23, 'A[A>T]A':9,  ...}},..}
	'''
	
	# Obtains all of the samples of interest from the input file
	out = open(log_file, 'a')
	sample_names = []
	samples = dict()
	for context in catalogue_files:
		with open (catalogue_files[context]) as f:
			first_line = f.readline().strip().split('\t')
		sample_names += first_line[1:]
		samples[context] = {}
		current_samples = first_line[1:]

		# Save the mutation counts for each sample for each nucleotide context
		with open (catalogue_files[context]) as f:
			next(f)
			for lines in f:
				line = lines.strip().split()
				nuc = line[0]
				if nuc == 'complex' or nuc == 'non-matching':
					continue

				sample_index = 1
				for sample in current_samples:
					mutCount = int(line[sample_index])
					if sample not in samples[context]:
						samples[context][sample] = {nuc:mutCount}
					else:
						samples[context][sample][nuc] = int(mutCount)
					sample_index += 1  

	sample_names = list(set(sample_names))
	print("Files successfully read and mutations collected. Mutation assignment starting now.", file=out)
	print("Files successfully read and mutations collected. Mutation assignment starting now.")
	out.close()
	return (sample_names, samples)	
	





def mut_tracker (sample_names, samples, reference_sample, nucleotide_context_files, chromosome_string_path, genome, chromosomes, bed, log_file):
	'''
	Returns a dictionary that contains the number of mutations allocated to each chromosome
		for a given nucleotide context for a given sample.
		
	Parameters:
				   sample_names  -> list of all samples
						samples  -> dictionary containing the mutation count for each mutation type
									for a given context.
			   reference_sample  -> uses the first sample in the list as a reference 
	   nucleotide_context_files  -> contains the chromosome proportion for each nucleotide for each context
		 chromosome_string_path  -> path to the chromosome reference files
						 genome  -> version of the genome desired as the reference
					chromosomes  -> list of chromosomes for the species of interest
							 bed -> flag that determines if the user has provided a BED file to simulate across
									a specific set of ranges

	Returns:
		mutation_tracker  -> a dictionary that contains the number of mutations allocated to each chromosome
							 for a given nucleotide context for each sample.

	Example return value:
		mutation_tracker = {'PDXXXX':{'A[A>C]A':{'X':2,'Y':1,'1':4,...},
									 {'A[A>T]A':{'X':2,'Y':1,...}, ...}
							'PDYYYY':{'A[A>C]A':{'X':1,'Y':3,'1':1,...},
									 {'A[A>T]A':{'X':3,'Y':2,...}, ...}}
	This function allocates the mutations based upon the size of the chromosome.
	
	'''
	out = open(log_file, 'a')
	mutation_tracker = {}

	### debugging ###
	# random_sample = random.sample(list(samples['384']),1)[0]
	# a = pd.DataFrame(0, index=list(samples['384'][random_sample].keys()), columns=chromosomes)
	# chrom_1_track_plus = 0
	# chrom_1_track_neg = 0
	# probs = [0.038016941,0.104438712,0.069723876,0.058404333,0.037029061,0.044805745,0.051422922,0.045144672,0.032261934,0.041033276,0.039722145,0.058713679,0.053670361,0.017297013,0.032746978,0.034216936,0.043138742,0.059986665,0.01449255,0.065151069,0.026191246,0.010340946,0.022050197]
	# probs = [0.05328559,0.079445171,0.084002805,0.06869521,0.066178886,0.062664235,0.059031866,0.054785456,0.050389775,0.042368533,0.046308117,0.04624282,0.046014254,0.033709769,0.03113531,0.028809649,0.027818706,0.027434474,0.026327851,0.019681034,0.020984617,0.012380337,0.012305537]
	
	#################
	for context in nucleotide_context_files:
		mutation_tracker[context] = {}
		sim = None
		mut_start = None
		mut_length = None
		if context == '6':
			sim = 1
		elif context == '24':
			sim = 10
		elif context == '96':
			sim = 2
		elif context == '1536':
			sim = 3
		elif context == '192' or context == '384':
			sim = 4 
		elif context == 'DINUC' or context == 'DBS':
			sim = 5
		elif context == 'INDEL' or context == 'ID':
			sim = 6
		elif context == '3072' or context == '6144':
			sim = 7
		elif context == 'DBS186' or context == '186':
			sim = 8
		elif context == 'ID415' or context == '415':
			sim = 9
		
		# Allocates mutations differently from INDEL simulation
		if sim != 6 and sim != 9:
			nuc_probs = {}

			# Opens and saves the context distribution for each nucleotide
			with open (nucleotide_context_files[context]) as f:
					next(f)
					for lines in f:
						line = lines.strip().split(',')
						nuc = line[0]
						line[1:] = list(map(float, line[1:]))
						nuc_probs[nuc] = line[1:]


			for sample in sample_names:
				random_sample = random.sample(list(samples[context]),1)[0]
				for nuc in samples[context][random_sample]:
					chrom_index = 0
					nuc_count = 0


					# Organizes the nucleotides and context dependent nucleotides
					# for allocation purposes
					if sim == 1:
						base_nuc = nuc[0]
					elif sim == 10:
						base_nuc = nuc[0:3]
					elif sim == 2: 
						base_nuc = nuc[0] + nuc[2] + nuc[6]
					elif sim == 4:
						base_nuc = nuc[0:3] + nuc[4] + nuc[8]
					elif sim == 3:
						base_nuc = nuc[0:2] + nuc[3] + nuc[7:]
					elif sim == 5:
						base_nuc = nuc[0:2]
					elif sim == 7:
						base_nuc = nuc[0:4] + nuc[5] + nuc[9:]
					elif sim == 8:
						base_nuc = nuc[0:4]


					# Pulls a distribution for randomly assigning left over mutations based upon the nucleotide distribution
					# across the genome
					probs = nuc_probs[base_nuc]
					if len(probs) != len(chromosomes):
						probs = probs[:len(chromosomes)]


					# Instantiate the mutation dictionary for the current nucleotide
					if sample not in mutation_tracker[context]:
						mutation_tracker[context][sample] = {nuc:{}}



					# Allocates mutations proportionaly based upon the context
					# distributions
					for chroms in chromosomes:							
						try:
							if sim == 5:
								try:
									nuc_probs[base_nuc][chrom_index]
								except:
									base_nuc = revcompl(base_nuc)
									nuc_probs[base_nuc][chrom_index]
							if sim == 8:
								try:
									nuc_probs[base_nuc][chrom_index]
								except:
									base_nuc = revbias(base_nuc[0]) + ":" + revcompl(base_nuc[2:])
									nuc_probs[base_nuc][chrom_index]								

							mutation_count = int(samples[context][sample][nuc]) * nuc_probs[base_nuc][chrom_index]

						except:
							mutation_count = 0
						if mutation_count - int(mutation_count) > 0.5:
							mutation_count = int(mutation_count) + 1
							nuc_count += mutation_count
						else:
							mutation_count = int(mutation_count)
							nuc_count += mutation_count
						if nuc not in mutation_tracker[context][sample]:
							mutation_tracker[context][sample][nuc] = {chroms:mutation_count}
						else:
							mutation_tracker[context][sample][nuc][chroms] = mutation_count
						chrom_index += 1

					# Ensures that the correct number of mutations have been assinged
					#if sample in samples[context].keys():
					if sample in samples[context]:
						if nuc_count != samples[context][sample][nuc]:
							while True:
								if nuc_count == samples[context][sample][nuc]:
									break
								else:
									# l = random.randint(0,len(chromosomes)-1)
									l = np.random.choice(len(chromosomes), p=probs)
									if nuc_probs[base_nuc][l] > 0:
										random_chromosome = chromosomes[l]

										if nuc_count < samples[context][sample][nuc]:
											mutation_tracker[context][sample][nuc][random_chromosome] += 1
											nuc_count += 1
										else:
											if mutation_tracker[context][sample][nuc][random_chromosome] != 0:
												mutation_tracker[context][sample][nuc][random_chromosome] -= 1
												nuc_count -= 1
					# for chrom in chromosomes:
					# 	a.at[nuc, chrom] += mutation_tracker[context][sample][nuc][chrom]

		# Allocates mutations for the INDEL simulation based upon the size
		# of each chromosome in relation to the overall size of the genome.
		else:
			if bed:
				with open (chromosome_string_path + "BED_" + genome + "_proportions.txt", 'rb') as probs:
					chromosomeProbs = pickle.load(probs)
			else:
				with open (chromosome_string_path + genome + "_proportions.txt", 'rb') as probs:
					chromosomeProbs = pickle.load(probs)

			for sample in sample_names:
				random_sample = random.sample(list(samples[context]),1)[0]
				for nuc in samples[context][random_sample]:
					chrom_index = 0;
					nuc_count = 0;
					if sample not in mutation_tracker[context]:
						mutation_tracker[context][sample] = {nuc:{}}

					# Allocates mutations based upong chromosome proportions
					for chroms in chromosomes:
						try:
							mutation_count = int(samples[context][sample][nuc]) * chromosomeProbs[chrom_index]
						except:
							mutation_count = 0
						if mutation_count - int(mutation_count) > 0.5:
							mutation_count = int(mutation_count) + 1
							nuc_count += mutation_count
						else:
							mutation_count = int(mutation_count)
							nuc_count += mutation_count
						if nuc not in mutation_tracker[context][sample]:
							mutation_tracker[context][sample][nuc] = {chroms:mutation_count}
						else:
							mutation_tracker[context][sample][nuc][chroms] = mutation_count
						chrom_index += 1

					# Ensures that the correct number of mutations have been assinged
					if sample in samples[context]:
						if nuc_count != samples[context][sample][nuc]:
							while True:
								if nuc_count == samples[context][sample][nuc]:
									break
								else:
									l = fastrand.pcg32bounded(21)
									random_chromosome = chromosomes[l]
									if nuc_count < samples[context][sample][nuc]:
										mutation_tracker[context][sample][nuc][random_chromosome] += 1
										nuc_count += 1
									else:
										if mutation_tracker[context][sample][nuc][random_chromosome] != 0:
											mutation_tracker[context][sample][nuc][random_chromosome] -= 1
											nuc_count -= 1


	print ("Mutations have been distributed. Starting simulation now...", file=out)
	print ("Mutations have been distributed. Starting simulation now...")
	out.close()
	return (mutation_tracker)
	
	
	
	
	
	
	
def simulator (sample_names, mutation_tracker, chromosome_string_path, tsb_ref, tsb_ref_rev, simulation_number, seed, cushion, output_path, updating, chromosomes, project, genome, bed, bed_file, contexts, overlap, project_path, seqInfo, log_file, spacing, noisePoisson, noiseAWGN, vcf):
	'''
	Simulates mutational signatures in human cancer in an unbiased fashion. The function
		requires a list of sample names, a dictionary of the number of mutations for each
		nucleotide context for each sample, and a dictionary with the mutations allocated
		proportionally to every chromosome. This function also requires that a local 
		copy of each chromosome be saved as a string within individual files ('1.txt', 
		'2.txt', '3.txt', etc). If TSB simulations are desired, the user must also save a
		local binary file for each chromosome that contains the transcriptional info (see
		blah.py for details on how to create the binary file for each chromosome). 
		
	Parameters:
		          sample_names  -> list of all samples 
						           (ex: sample_names = ['PDXXXX', 'PDYYYY', ...])
		
					   samples  -> dicationary with mutation counts for each mutation type for each
								   sample.
								   (ex: samples  -> {'PDXXXX':{'A[A>C]A':35, 'A[A>T]A':12, ...}
													 'PDYYYY':{'A[A>C]A':23, 'A[A>T]A':9,  ...}})
			  mutation_tracker  -> {'PDXXXX':{'A[A>C]A':{'X':2,'Y':1,'1':4,...},
											 {'A[A>T]A':{'X':2,'Y':1,...}, ...}
									'PDYYYY':{'A[A>C]A':{'X':1,'Y':3,'1':1,...},
											 {'A[A>T]A':{'X':3,'Y':2,...}, ...}}
		chromosome_string_path  -> path to the chromosome reference files
		               tsb_ref  -> dictionary that allows switching from binary code to biologically relevant strings
		           tsb_ref_rev  -> dictionary that allows switching from biologically relevant strings back to binary values
			 simulation_number  -> desired simulation number
				   output_path  -> output path for the simulations
					  updating  -> single value to determine whether updating should occur. 
				   chromosomes  -> list of chromosomes for the given genome
					   project  -> unique name for the given samples
						genome  -> reference genome used for simulation
						   bed  -> flag to determine if a BED file with a specific set of ranges was provided
					  bed_file  -> if specified by the user, the BED file with the given set of ranges. Else,
								   it will be equal to 'None'
					  contexts  -> desired nucleotide contexts for simulation
					     exome  -> flag that simulates based upon the exome
					   overlap  -> flag that allows SNV mutations and DBS mutations to overlap. By default, they will not overlap.

	Returns:
		None 

	Output: 
		Writes the output to a single vcf file per folder for each desired context.
		See https://samtools.github.io/hts-specs/VCFv4.2.pdf for an example vcf format.

	'''
	# Saves the chromosome as a string in memory and creates a list if a BED file was supplied
	# tracker = SummaryTracker()

	left_over_mutations = {}
	
	dinuc_non_tsb = ['AC', 'AT', 'CA', 'CG', 'GC', 'TA']

	# Set seed
	try:
		fastrand.pcg32_seed(seed)
	except:
		seed = seed % 100
		fastrand.pcg32_seed(seed)
		
	if seqInfo:
		seqOut_path = project_path + "output/vcf_files/simulations/"


	ranges = {}
	file_context = "_".join(contexts)
	for chrom in chromosomes:
		if bed:
			chrom_range = bed_ranges (chrom, bed_file, cushion)
		with open (chromosome_string_path + chrom + ".txt", "rb") as f:
			initial_seq = f.read().strip()
			if updating:
				initial_seq = list(initial_seq)
				
		# Only for TSB simulations, opens the transcriptional info strings:
		chrom_bias = {'TU':[],'B':[],'N':[]}
		chrom_bias_lengths = {'TU':[],'B':[],'N':[]}
		if '192' in contexts or '3072' in contexts or '384' in contexts or '6144' in contexts or 'DBS186' in contexts or '24' in contexts or 'ID415' in contexts:
			chromosome_string_path, ref_dir = matRef.reference_paths(genome)
			with open (ref_dir + '/references/chromosomes/tsb_BED/' + genome + '/' + chrom + "_BED_TSB.txt") as f:
				next(f)
				for lines in f:
					line = lines.strip().split()
					start = int(line[1])
					end = int(line[2])
					range_bias = int(line[3])
					if range_bias == 0:
						chrom_bias['N'].append([start, end])
						if chrom_bias_lengths['N'] == []:
							chrom_bias_lengths['N'].append(end-start)
						else:
							chrom_bias_lengths['N'].append(chrom_bias_lengths['N'][-1] + (end-start))
					elif range_bias == 1:
						chrom_bias['TU'].append([start, end])
						if chrom_bias_lengths['TU'] == []:
							chrom_bias_lengths['TU'].append(end-start)
						else:
							chrom_bias_lengths['TU'].append(chrom_bias_lengths['TU'][-1] + (end-start))
					elif range_bias == 2:
						chrom_bias['TU'].append([start, end])
						if chrom_bias_lengths['TU'] == []:
							chrom_bias_lengths['TU'].append(end-start)
						else:
							chrom_bias_lengths['TU'].append(chrom_bias_lengths['TU'][-1] + (end-start))
					elif range_bias == 3:
						chrom_bias['B'].append([start, end])
						if chrom_bias_lengths['B'] == []:
							chrom_bias_lengths['B'].append(end-start)
						else:
							chrom_bias_lengths['B'].append(chrom_bias_lengths['B'][-1] + (end-start))

		for sample in sample_names:
			simulations = simulation_number
			sample_path = output_path
			if vcf:
				sample_path = sample_path + sample + "/"

			while(simulations > 0):
				
				# Saves an unaltered chromosome, so that future updating of mutations
				# does not affect additional simulations.
				sequence = initial_seq
				if bed:
					location_range = len(chrom_range)
					if location_range == 0:
						break
				else:
					location_range = len(sequence)
				recorded_positions = set()

				# Creates the output path if it does not already exist.  
				if not os.path.exists(output_path):
					os.makedirs(output_path)  

				if vcf:
					outputFile = ''.join([sample_path,sample,"_",str(simulations),"_",chrom,".vcf"])
				else:
					outputFile = ''.join([sample_path,str(simulations),"_",chrom,".maf"])

				with open(outputFile, "a", 10000000) as out_vcf:
					for context in contexts:
						if seqInfo:
							outSeq = open(seqOut_path + context + "/" + sample + "_" + chrom + "_seqinfo_" + str(simulations) + ".txt", "w", 10000000) 

						sim = None
						mut_start = None
						mut_length = None
						if context == '6':
							sim = 1
							mut_start = 0
							mut_save = 0
						elif context == '24':
							sim = 10
							mut_start = 0
							mut_save = 2
						elif context == '96':
							sim = 2
							mut_start = 1
							mut_save = 2
						elif context == '1536':
							sim = 3
							mut_start = 2
							mut_save = 3
						elif context == '192' or context == '384':
							sim = 4 
							mut_start = 1
							mut_save = 4
						elif context == 'DINUC' or context == 'DBS':
							sim = 5
							mut_start = 0
							mut_save = 0
						elif context == 'INDEL' or context == 'ID':
							sim = 6
							mut_save = 0
							mut_start = 0
						elif context == '3072' or context == '6144':
							sim = 7
							mut_save = 5
							mut_start = 2
						elif context == 'DBS186':
							sim = 8
							mut_start = 1
							mut_save = 5	
						elif context == 'ID415':
							sim = 9
							mut_start = 0
							mut_save = 0						
						
						mutationsCount = {}
						random_sample = random.sample(list(mutation_tracker[context]),1)[0]
					
						for nuc in mutation_tracker[context][random_sample]:
							mutationsCount[nuc] = mutation_tracker[context][sample][nuc][chrom]
						if sample in left_over_mutations and simulations in left_over_mutations[sample]: 

							if any(left_over_mutations[sample][simulations]):
								for nuc in left_over_mutations[sample][simulations][context]:
									try:
										mutationsCount[nuc] += left_over_mutations[sample][simulations][context][nuc]
									except:
										mutationsCount[nuc] = left_over_mutations[sample][simulations][context][nuc]
									left_over_mutations[sample][simulations][context][nuc] = {}

						# Add in noise:
						if noisePoisson or noiseAWGN != 0:
							mutationsCount = noise(mutationsCount, noisePoisson, noiseAWGN)



						initial_nuc_keys = list(mutationsCount.keys())
						for nuc in initial_nuc_keys:
							if mutationsCount[nuc] == 0:
								del mutationsCount[nuc]

						nuc_keys = list(mutationsCount.keys())
						base_keys = []


						# Simulations for INDEL context
						if sim == 6:
							u = 0
							indel_lengths = []
							repeat_lengths = []
							indel_types = {}
							indel_types_O = {}
							indel_types_M = {}

							# Organizes data structures to keep track of the desired INDEL mutations
							for indels in mutationsCount:
								indel = indels.split(':')
								if sim == 9:
									bias = indel[0]
									indel = indel[1:]
								if int(indel[0]) == 5 and int(indel[3]) == 5 and indel[2] == 'M':
									indel_lengths.append(6)
								else:
									indel_lengths.append(int(indel[0]))

								repeat_lengths.append(int(indel[3]))
								save_key = ''
								if sim == 9:
									save_key = bias
								# Saves the Insertion mutations with 0 repeats in a separate data structure
								if indel[3] == '0' and indel[1] == 'Ins':
									if (indel[0]+indel[3] + indel[2] + save_key) not in indel_types_O:
										indel_types_O[(indel[0]+indel[3] + indel[2] + save_key)] = [indel[1]]

								# Saves the Insertion mutations for microhomologies separately
								elif indel[1] == 'Ins' and indel[2] == 'M':
									if (indel[0]+indel[3] + indel[2]+ save_key) not in indel_types_M:
										indel_types_M[(indel[0]+indel[3] + indel[2] + save_key)] = [indel[1]]

								# Saves all other INDEL mutations in a single data structure
								else:
									if (indel[0]+indel[3] + indel[2]+ save_key) not in indel_types:
										indel_types[(indel[0]+indel[3] + indel[2] + save_key)] = [indel[1]]
									else:
										indel_types[(indel[0]+indel[3] + indel[2] + save_key)].append(indel[1])

							# Repeats simulation until all mutations are assigned
							while (any(mutationsCount) == True):
								while (any(indel_types)) == True:
									u += 1
									if u > 1000000:
										# logging.info(sample + " " + mutationsCount)
										print (sample, mutationsCount, chrom)
										indel_types = {}
										u = 0

									break_flag = False
									random_number = fastrand.pcg32bounded(location_range)
									if bed:
										random_number = chrom_range[random_number]

									if random_number in recorded_positions and not overlap:
										continue
									u += 1
									stop_flag = False

									# Pulls out the largest desired INDEL
									for i in range (max(indel_lengths), 0, -1):
										micro = False
										inDel = ''
										try:
											for r in range (random_number,i+random_number,1):
												if r > len(sequence):
													break
												inDel += tsb_ref[sequence[r]][1]
										except:
											break


										# Ensures that all bases are known in the potential mutation
										if 'N' not in inDel:
											repeat_count = 0
											repeat_count_ins = 0

											# Counts the number of repeats of the INDEL in the forward direction
											for k in range (1, 6, 1):
												seq = ''
												try:
													for r in range(random_number+(k*i), random_number+((k+1)*i), 1):
														if r > len(sequence):
															break
														seq += tsb_ref[sequence[r]][1]
												except:
													break   

												if seq == inDel:
													repeat_count += 1
												else:
													if repeat_count == 0 and len(inDel) > 1:
														for p in range(1, len(seq), 1):
															if seq[:p] == inDel[:p]:
																micro=True
																break
													break


											# Counts the number of repeats of the INDEL in the reverse direction
											for l in range (1, 6, 1):
												seq = ''
												try:
													for r in range(random_number-(l*i),(random_number-(l*i))+i, 1):
														if r > len(sequence):
															break
														seq += tsb_ref[sequence[r]][1]
												except:
													break
												if seq == inDel:
													repeat_count += 1
													micro=False
												else:
													if repeat_count == 0 and len(inDel) > 1:
														for p in range(1, len(seq), 1):
															if seq[-p:] == inDel[-p:]:
																micro=True
																break
													break

											# Organizes a naming structure for the current INDEL chosen on the chromosome
											mainType = str(i) + str(repeat_count+repeat_count_ins)
											mainType_ins = str(i) + str(repeat_count+1+repeat_count_ins) 


											# Assigns the subtype category for the INDEL based on its context
											strand = "+"
											if i != 1:
												subType = 'R'
												if bool(re.match("^[CT]*$", inDel)) or bool(re.match("^[GA]*$", inDel)):
													if bool(re.match("^[GA]*$", inDel)):
														strand = "-"
												else:
													strand = "Q"
											else:
												if tsb_ref[sequence[random_number]][1] == 'A' or tsb_ref[sequence[random_number]][1] == 'T':
													subType = 'T'
													if tsb_ref[sequence[random_number]][1] == 'A':
														strand = "-"

												else:
													subType = 'C'
													if tsb_ref[sequence[random_number]][1] == 'G':
														strand = "-"

											# Saves a type for a deleltion and an insertion option
											mainType += subType
											mainType_ins += subType

											if micro:
												mainType = None
												if mainType_ins not in indel_types:
													mainType_ins = None

											# Checks to see if this randomly chosen INDEL is a desired deletion mutation
											if mainType in indel_types:
												if indel_types[mainType][0] == 'Del':
													if not overlap:
														for p in range(random_number, random_number+i, 1):
															if p in recorded_positions:
																stop_flag = True
													if stop_flag:
														break

													# Reassigns the original naming convention for the INDEL and writes it to the output file
													complete_indel = mainType[0] + ':Del:' + subType + ':' + mainType[1]
													seq_final = ''
													for r in range (random_number-1,i+random_number,1):
														seq_final += tsb_ref[sequence[r]][1]

													if vcf:
														print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
													else:
														print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",seq_final,".",tsb_ref[sequence[random_number-1]][1],".\t.",sample+"_"+str(simulations),complete_indel]), file=out_vcf)


													if seqInfo:
														print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)

													if not overlap:
														for z in range (random_number-1, random_number+i+1,1):
															recorded_positions.add(z)

													# Updates the chromosome with the given mutation if desired
													if updating:
														sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Del')
														location_range = len(sequence)

													# Remove one mutation count for the current INDEL and update all relevant data structures
													mutationsCount[complete_indel] -= 1
													u = 0
													if mutationsCount[complete_indel] == 0:
														del mutationsCount[complete_indel]
														indel_lengths.remove(int(mainType[0]))
														repeat_lengths.remove(repeat_count+repeat_count_ins)
														if len(indel_types[mainType]) > 1:
															del indel_types[mainType][0]
														else:
															del indel_types[mainType]
													break

											# Checks to see this randomly chosen INDEL is a desired insertion mutation
											elif mainType_ins in indel_types:
												if indel_types[mainType_ins][0] == 'Ins':
													if not overlap:
														for p in range(random_number, random_number+i, 1):
															if p in recorded_positions:
																stop_flag = True
													if stop_flag:
														break

													#Reassigns the original naming convention for the INDEL and writes it to the output file
													complete_indel = mainType_ins[0] + ':Ins:' + subType + ':' + mainType_ins[1]
													potential_sequence = ''
													for r in range(random_number, random_number+i,1):
														potential_sequence += tsb_ref[sequence[r]][1]
													seq_final = ''
													for r in range (random_number-1, i+random_number,1):
														seq_final += tsb_ref[sequence[r]][1]
													
													if vcf:
														print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
													else:
														print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",tsb_ref[sequence[random_number-1]][1],".",seq_final,".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)

													if seqInfo:
														print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)

													if not overlap:
														for z in range (random_number-1, random_number+i+1,1):
															recorded_positions.add(z)

													# Updates the chromosome with the given mutation if desired
													if updating:
														sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Ins')
														location_range = len(sequence)


													# Remove one mutation count for the current INDEL and update all relevant data structures
													try:
														mutationsCount[complete_indel] -= 1
														# print(chrom, mutationsCount)
														u = 0
													except:
														sys.stderr.write(mainType_ins, complete_indel, indel_types, mutationsCount)
														print(mainType_ins, complete_indel, indel_types, mutationsCount)
													if mutationsCount[complete_indel] == 0:
														del mutationsCount[complete_indel]
														indel_lengths.remove(int(mainType_ins[0]))
														repeat_lengths.remove(repeat_count+1+repeat_count_ins)
														if len(indel_types[mainType_ins]) > 1:
															del indel_types[mainType_ins][0]
														else:
															del indel_types[mainType_ins]
													break


											# Simulates microhomology deletions
											else:
												if repeat_count != 0 or i == 1:
													break
												max_repeat_length = max(repeat_lengths)+1
												if max_repeat_length > i:
													max_repeat_length = i


												# Counts homology in the forward direction
												homology_size1 = 0
												for k in range (1, 6, 1):
													seq = ''
													try:
														for r in range (random_number+i, random_number+k+i, 1):
															if r > len(sequence):
																break
															seq += tsb_ref[sequence[r]][1]
													except:
														break
													if seq == inDel[:k]:
														homology_size1 += 1
													else:
														break

												# Counts homology in the reverse direction
												homology_size2  = 0
												for l in range (1, 6, 1):
													seq = ''
													try:
														for r in range (random_number-l, random_number, 1):
															if r > len(sequence):
																break
															seq += tsb_ref[sequence[r]][1]
													except:
														break
													if seq == inDel[-l:]:
														homology_size2 += 1
													else:
														break

												# Assigns a new naming convetion for the current INDEL
												subType = 'M'
												if i > 5 and (homology_size1 >= 5 or homology_size2 >= 5):
													mainType1 = str(i-1) + str(homology_size1) + subType
													mainType2 = str(i-1) + str(homology_size2) + subType
												else:
													mainType1 = str(i) + str(homology_size1) + subType
													mainType2 = str(i) + str(homology_size2) + subType
												complete_indel = None
												
												# Checks to see if the forward homology is desired and that the reverse
												# homology equals 0. 
												if mainType1 in indel_types and homology_size2 == 0:
														if not overlap:
															for p in range(random_number, random_number+i, 1):
																if p in recorded_positions:
																	stop_flag = True
														if stop_flag:
															break

														# Reassigns the original naming convention for the INDEL and writes it to the output file
														complete_indel = mainType1[0] + ':Del:' + subType + ':' + mainType1[1]
														seq_final = ''
														for r in range(random_number-1,i+random_number,1):
															seq_final += tsb_ref[sequence[r]][1]
														if bool(re.match("^[CT]*$", seq_final)) or bool(re.match("^[GA]*$", seq_final)):
															if bool(re.match("^[GA]*$", seq_final)):
																strand = "-"
														else:
															strand = "Q"
														if vcf:
															print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
														else:
															print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",seq_final,".",tsb_ref[sequence[random_number-1]][1],".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
														
														if seqInfo:
															print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
														
														if not overlap:
															for z in range (random_number-1, random_number+i+1,1):
																recorded_positions.add(z)

														# Updates the chromosome with the current INDEL if desired
														if updating:
															sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Del')
												 
														# Remove one mutation count for the current INDEL and update all relevant data structures       
														mutationsCount[complete_indel] -= 1
														u = 0
														if mutationsCount[complete_indel] == 0:
															del mutationsCount[complete_indel]
															indel_lengths.remove(i)
															repeat_lengths.remove(homology_size1)
															if len(indel_types[mainType1]) > 1:
																del indel_types[mainType1][0]
															else:
																del indel_types[mainType1]
														break


												# Checks to see if the reverse homology is desired and that the forward
												# homology equals 0.
												elif mainType2 in indel_types and homology_size1 == 0:
													seq = ''
													try:
														for r in range (random_number+i, random_number+(2*i),1):
															if r > len(sequence):
																break
															seq += tsb_ref[sequence[r]][1]
													except:
														break

													seq2 = ''
													try:
														for r in range(random_number, random_number+i,1):
															if r > len(sequence):
																break
															seq2 += tsb_ref[sequence[r]][1]
													except:
														break

													seq3 = ''
													try:
														for r in range(random_number-i, random_number, 1):
															if r > len(sequence):
																break
															seq3 += tsb_ref[sequence[r]][1]
													except:
														break

													if indel_types[mainType2][0] == 'Del' and seq != seq2 and seq3 != seq2:
														if not overlap:
															for p in range(random_number, random_number+i, 1):
																if p in recorded_positions:
																	stop_flag = True
														if stop_flag:
															break

														# Reassigns the original naming convention for the INDEL and writes it to the output file
														complete_indel = mainType2[0] + ':Del:' + subType + ':' + mainType2[1]
														
														seq_final = ''
														for r in range (random_number-1,i+random_number,1):
															seq_final += tsb_ref[sequence[r]][1]
														if bool(re.match("^[CT]*$", seq_final)) or bool(re.match("^[GA]*$", seq_final)):
															if bool(re.match("^[GA]*$", seq_final)):
																strand = "-"
														else:
															strand = "Q"

														if vcf:
															print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
														else:
															print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",seq_final,".",tsb_ref[sequence[random_number-1]][1],".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
										
														if seqInfo:
															print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
														
														if not overlap:
															for z in range (random_number-1, random_number+i+1,1):
																recorded_positions.add(z)

														# Updates the chromosome with the current INDEL if desired
														if updating:
															sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Del')
															location_range = len(sequence)

														# Remove one mutation count for the current INDEL and update all relevant data structures
														mutationsCount[complete_indel] -= 1
														u = 0
														if mutationsCount[complete_indel] == 0:
															del mutationsCount[complete_indel]
															indel_lengths.remove(i)
															repeat_lengths.remove(homology_size2)
															if len(indel_types[mainType2]) > 1:
																del indel_types[mainType2][0]
															else:
																del indel_types[mainType2]
														break


								# Simuales all micro-homology insertion mutations
								for indels_M in indel_types_M:
									if u > 1000000:
										print (sample, mutationsCount, chrom)
										indel_types_M = {}
										u = 0
									# Randomly chooses a location on the current chromosome
									random_number = fastrand.pcg32bounded(location_range)
									u += 1

									# Assigns the original naming convention
									complete_indel = indels_M[0] + ':Ins:' + indels_M[2] + ':' + indels_M[1]
									
									# Pulls the sequence of bases out for the current insertion homology
									# equal to the length of the homology
									M_length = int(indels_M[0])
									if int(indels_M[0]) == 5 and int(indels_M[1]) == 5:
										M_length = 6

									potential_sequence = ''
									for r in range (random_number, random_number + int(indels_M[1]), 1):
										if r > len(sequence):
											break
										potential_sequence += tsb_ref[sequence[r]][1]

									if not bool(re.match("^[CT]*$", inDel)) and not bool(re.match("^[GA]*$", inDel)):
										bias = 'Q'
										strand = "Q"
									else:
										if bool(re.match("^[GA]*$", inDel)):
											strand = "-"
										bias = tsb_ref[sequence[random_number]][0]

									# Esnures that the region is known
									if 'N' not in potential_sequence:
										
										# Saves the potential reverese homology sequence for reference. 
										reverse_homology = ''
										for r in range (random_number-int(indels_M[1]), random_number, 1):
											if r > len(sequence):
												break
											reverse_homology += tsb_ref[sequence[r]][1]
										remaining_sequence = ''

										# Adds random bases to the end of the micro-homology, ensuring
										# that the added bases don't exceed the desired homology length
										for m in range (0, M_length-int(indels_M[1])-1, 1):
											while len(remaining_sequence) != 1:
												new_base = random_base(0,3)
												if new_base != tsb_ref[sequence[random_number+int(indels_M[1])]][1]:
													remaining_sequence += new_base
											potential_sequence += new_base#random_base(0,3)


										# Adds random bases until the desired insertion length is met,
										# without introducing reverse homology
										while len(potential_sequence) != M_length:
											last_base = random_base(0,3)
											if last_base != reverse_homology[-1]:
												potential_sequence += last_base

										# Prints the insertion micro-hommology if the insertion length is correct and the homology matches are correct
										seq = ''
										try:
											for r in range(random_number-int(indels_M[1]),random_number,1):
												if r > len(sequence):
													break
												seq += tsb_ref[sequence[r]][1]
											seq2 = ''
											for r in range (random_number+M_length+int(indels_M[1]),random_number+M_length+(2*int(indels_M[1])), 1):
												if r > len(sequence):
													break
												seq2 += tsb_ref[sequence[r]][1]
											seq3 = ''
											for r in range(random_number,random_number+int(indels_M[1])+1,1):
												if r > len(sequence):
													break
												seq3 += tsb_ref[sequence[r]][1]
										except:
											break

										if seq != potential_sequence[-int(indels_M[1]):] and seq2 != potential_sequence[:int(indels_M[1])] and seq3 != potential_sequence[:int(indels_M[1])+1]:
											if not overlap:
												for p in range(random_number, random_number+M_length+1, 1):
													if p in recorded_positions:
														stop_flag = True
											if stop_flag:
												break

											if vcf:
												print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", tsb_ref[sequence[random_number-1]][1],"\t",tsb_ref[sequence[random_number-1]][1]+potential_sequence,"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
											else:
												print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",tsb_ref[sequence[random_number-1]][1],".",tsb_ref[sequence[random_number-1]][1]+potential_sequence,".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
											
											if seqInfo:
												print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
											
											if not overlap:
												for z in range (random_number-1, random_number+M_length+1,1):
													recorded_positions.add(z)

											# Remove one mutation count for the current INDEL and update all relevant data structures
											mutationsCount[complete_indel] -= 1
											u = 0
											if mutationsCount[complete_indel] == 0:
												del mutationsCount[complete_indel]
												indel_lengths.remove(M_length)
												repeat_lengths.remove(int(indels_M[1]))
												del indel_types_M[indels_M]
											break

											# Updates the chromosome with the current INDEL if desired
											if updating:
												sequence = update_chromosome(sequence, random_number, sequence[random_number:int(indels_M[0])+random_number], 'Ins')
												location_range = len(sequence)


								# Simulates the insertion INDELs that have 0 repeats
								if not any(indel_types_O):
									mutationsCount = {}
								while (any(indel_types_O) == True):
									if u > 1000000:
										# logging.info(sample + " " + mutationsCount)
										print (sample, mutationsCount, chrom)
										indel_types_O = {}
										mutationsCount = {}
										u = 0
									# Randomly chooses a location on the current chromosome
									random_number = fastrand.pcg32bounded(location_range)
									if random_number in recorded_positions and not overlap:
										continue	
									u += 1

									# Assigns a subtype for the chosen position on the chromosome
									strand = "+"
									if tsb_ref[sequence[random_number]][1] == 'T' or tsb_ref[sequence[random_number]][1] == 'A':
										subType = 'T'
									elif tsb_ref[sequence[random_number]][1] == 'G' or tsb_ref[sequence[random_number]][1] == 'C':
										subType = 'C'
									else:
										break

									stop_flag = False
									for indels_O in indel_types_O:

										# For INDELs of length 1, if the subType for the mutation does not match
										# the desired INDEL, break and find a new position
										if indels_O[0] == '1':
											if indels_O[2] != subType:
												break

										# Assigns the original naming convention
										complete_indel = indels_O[0] + ':Ins:' + indels_O[2] + ':' + indels_O[1]

										# Randomly assigns a base for insertions of length 1 based upon
										# the subType at the chosen position in the chromosome.
										potential_sequence = ''
										if int(indels_O[0]) == 1:
											if subType == 'T':
												potential_sequence += random_base(0,1)
												if tsb_ref[sequence[random_number]][1] == 'A':
													strand = "-"
											else:
												potential_sequence += random_base(2,3)
												if tsb_ref[sequence[random_number]][1] == 'G':
													strand = "-"

										# Randomly chooses bases until the insertion length equals the desired
										# INDEL length
										else:
											for m in range (0, int(indels_O[0]), 1):
												potential_sequence += random_base(0,3)


										# Ensures that the bases are known and that there are no repeats around the insertion
										try:
											seq = ''
											for r in range (random_number-int(indels_O[0]),random_number,1):
												if r > len(sequence):
													break
												seq += tsb_ref[sequence[r]][1]
											seq2 = ''
											for r in range (random_number,random_number+int(indels_O[0]),1):
												if r > len(sequence):
													break
												seq2 += tsb_ref[sequence[r]][1]
										except:
											break

										if "N" not in seq:
											if seq2 != potential_sequence and tsb_ref[sequence[random_number-int(indels_O[0])]][1] != potential_sequence:
												if not overlap:
													for p in range(random_number, random_number+int(indels_O[0])+1, 1):
														if p in recorded_positions:
															stop_flag = True
												if stop_flag:
													break
													
												if vcf:
													print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", tsb_ref[sequence[random_number-1]][1],"\t",tsb_ref[sequence[random_number-1]][1]+potential_sequence,"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
												else:
													print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",tsb_ref[sequence[random_number-1]][1],".",tsb_ref[sequence[random_number-1]][1]+potential_sequence,".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
												if seqInfo:
													print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
												if not overlap:
													for z in range (random_number-1, random_number+int(indels_O[0])+1,1):
														recorded_positions.add(z)

												# Remove one mutation count for the current INDEL and update all relevant data structures
												mutationsCount[complete_indel] -= 1
												# print(chrom, mutationsCount)
												u = 0
												if mutationsCount[complete_indel] == 0:
													del mutationsCount[complete_indel]
													indel_lengths.remove(int(indels_O[0]))
													repeat_lengths.remove(int(indels_O[1]))
													del indel_types_O[indels_O]
												break

												# Updates the chromosome with the current INDEL if desired
												if updating:
													sequence = update_chromosome(sequence, random_number, sequence[random_number:int(indels_O[0])+random_number], 'Ins')
													location_range = len(sequence)
										break
						

						# Simulations for INDEL with TSB
						elif sim == 9:
							tsb = ['TU','B','N', 'Q']
							mutationsCountTSB = {'TU':{},'B':{},'N':{},'Q':{}}
							indel_lengths = {'TU':[],'B':[],'N':[], 'Q':[]}
							repeat_lengths = {'TU':[],'B':[],'N':[], 'Q':[]}
							indel_types = {'TU':{},'B':{},'N':{}, 'Q':{}}
							indel_types_O = {'TU':{},'B':{},'N':{}, 'Q':{}}

							# Organizes data structures to keep track of the desired INDEL mutations
							for indels in mutationsCount:
								indel = indels.split(':')
								bias = indel[0]
								if bias == 'T' or bias == 'U':
									refbias = 'TU'
								else:
									refbias = bias
								indel = indel[1:]
								mutationsCountTSB[refbias][indels] = mutationsCount[indels]
								if int(indel[0]) == 5 and int(indel[3]) == 5 and indel[2] == 'M':
									indel_lengths[refbias].append(6)
								else:
									indel_lengths[refbias].append(int(indel[0]))
								repeat_lengths[refbias].append(int(indel[3]))

								# Saves the Insertion mutations with 0 repeats in a separate data structure
								if indel[3] == '0' and indel[1] == 'Ins':
									if (indel[0]+indel[3] + indel[2] + bias) not in indel_types_O[refbias]:
										indel_types_O[refbias][(indel[0]+indel[3] + indel[2] + bias)] = [indel[1]]


								# Saves all other INDEL mutations in a single data structure
								else:
									if (indel[0]+indel[3] + indel[2]+ bias) not in indel_types[refbias]:
										indel_types[refbias][(indel[0]+indel[3] + indel[2] + bias)] = [indel[1]]
									else:
										indel_types[refbias][(indel[0]+indel[3] + indel[2] + bias)].append(indel[1])

							# Repeats simulation until all mutations are assigned
							for tsb_type in tsb:
								u = 0
								while (any(mutationsCountTSB[tsb_type]) == True):
									while (any(indel_types[tsb_type])) == True:
										if u > 1000000:
											# print(indel_types[tsb_type])
											indel_types[tsb_type] = {}
											u = 0

										break_flag = False
										u += 1
										if not bed:
											if tsb_type != 'Q':
												location_range = chrom_bias_lengths[tsb_type][-1]
												random_range = fastrand.pcg32bounded(location_range)
												specific_range = bisect.bisect_left(chrom_bias_lengths[tsb_type], random_range)
												random_number = (chrom_bias_lengths[tsb_type][specific_range] - random_range) + chrom_bias[tsb_type][specific_range][0]
											else:
												random_number = fastrand.pcg32bounded(location_range)

										else:
											if tsb_type != 'Q':
												location_range = len(chrom_range)
												random_range = fastrand.pcg32bounded(location_range)
												random_number = chrom_range[random_range]
											else:
												random_number = fastrand.pcg32bounded(location_range)
												random_number = chrom_range[random_number]

										if random_number in recorded_positions and not overlap:
											continue

										stop_flag = False

										
										# Pulls out the largest desired INDEL
										for i in range (max(indel_lengths[tsb_type]), 0, -1):
											strand = "+"
											micro = False
											bias = tsb_ref[sequence[random_number]][0]
											inDel = ''
											try:
												for r in range (random_number,i+random_number,1):
													if r > len(sequence):
														break
													inDel += tsb_ref[sequence[r]][1]
											except:
												break

											if bias != 'N':
												if not bool(re.match("^[CT]*$", inDel)) and not bool(re.match("^[GA]*$", inDel)):
													bias = 'Q'
													if tsb_type != bias:
														continue
												else:
													if bool(re.match("^[GA]*$", inDel)):
														bias = revbias(bias)

											if bool(re.match("^[GA]*$", inDel)) or bool(re.match("^[CT]*$", inDel)):
												if bool(re.match("^[GA]*$", inDel)):
													strand = "-"
											else:
												strand = "Q"
											# Ensures that all bases are known in the potential mutation
											if 'N' not in inDel:
												repeat_count = 0
												repeat_count_ins = 0


												# Counts the number of repeats of the INDEL in the forward direction
												for k in range (1, 6, 1):
													seq = ''
													try:
														for r in range(random_number+(k*i), random_number+((k+1)*i), 1):
															if r > len(sequence):
																break
															seq += tsb_ref[sequence[r]][1]
													except:
														break   

													if seq == inDel:
														repeat_count += 1
													else:
														if repeat_count == 0 and len(inDel) > 1:
															for p in range(1, len(seq), 1):
																if seq[:p] == inDel[:p]:
																	micro=True
																	break
														break


												# Counts the number of repeats of the INDEL in the reverse direction
												for l in range (1, 6, 1):
													seq = ''
													try:
														for r in range(random_number-(l*i),(random_number-(l*i))+i, 1):
															if r > len(sequence):
																break
															seq += tsb_ref[sequence[r]][1]
													except:
														break
													if seq == inDel:
														repeat_count += 1
														micro=False
													else:
														if repeat_count == 0 and len(inDel) > 1:
															for p in range(1, len(seq), 1):
																if seq[-p:] == inDel[-p:]:
																	micro=True
																	break
														break

												# Organizes a naming structure for the current INDEL chosen on the chromosome
												mainType = str(i) + str(repeat_count+repeat_count_ins)
												mainType_ins = str(i) + str(repeat_count+1+repeat_count_ins) 


												# Assigns the subtype category for the INDEL based on its context
												if i != 1:
													subType = 'R'
												else:
													if tsb_ref[sequence[random_number]][1] == 'A' or tsb_ref[sequence[random_number]][1] == 'T':
														subType = 'T'
													else:
														subType = 'C'
											
												# Saves a type for a deleltion and an insertion option
												mainType += subType
												mainType_ins += subType

												# Saves bias if simulating TSB
												mainType += bias
												mainType_ins += bias
												# print(mainType_ins, indel_types[tsb_type])
												if micro:
													mainType = None
													if mainType_ins not in indel_types:
														mainType_ins = None

												# Checks to see if this randomly chosen INDEL is a desired deletion mutation
												if mainType in indel_types[tsb_type]:
													if indel_types[tsb_type][mainType][0] == 'Del':
														if not overlap:
															for p in range(random_number, random_number+i, 1):
																if p in recorded_positions:
																	stop_flag = True
														if stop_flag:
															break

														# Reassigns the original naming convention for the INDEL and writes it to the output file
														complete_indel = mainType[-1] + ':' +mainType[0] + ':Del:' + subType + ':' + mainType[1]
														seq_final = ''
														for r in range (random_number-1,i+random_number,1):
															seq_final += tsb_ref[sequence[r]][1]

														if vcf:
															print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
														else:
															print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",seq_final,".",tsb_ref[sequence[random_number-1]][1],".\t.",sample+"_"+str(simulations),complete_indel]), file=out_vcf)

														if seqInfo:
															print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)

														if not overlap:
															for z in range (random_number-1, random_number+i+1,1):
																recorded_positions.add(z)

														# Updates the chromosome with the given mutation if desired
														if updating:
															sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Del')
															location_range = len(sequence)

														# Remove one mutation count for the current INDEL and update all relevant data structures
														mutationsCountTSB[tsb_type][complete_indel] -= 1
														
														u = 0
														if mutationsCountTSB[tsb_type][complete_indel] == 0:
															del mutationsCountTSB[tsb_type][complete_indel]
															indel_lengths[tsb_type].remove(int(mainType[0]))
															repeat_lengths[tsb_type].remove(repeat_count+repeat_count_ins)
															if len(indel_types[tsb_type][mainType]) > 1:
																del indel_types[tsb_type][mainType][0]
															else:
																del indel_types[tsb_type][mainType]
														# print(chrom, tsb_type, mutationsCountTSB[tsb_type])
														break

												# Checks to see this randomly chosen INDEL is a desired insertion mutation
												elif mainType_ins in indel_types[tsb_type]:
													if indel_types[tsb_type][mainType_ins][0] == 'Ins':
														if not overlap:
															for p in range(random_number, random_number+i, 1):
																if p in recorded_positions:
																	stop_flag = True
														if stop_flag:
															break

														#Reassigns the original naming convention for the INDEL and writes it to the output file
														complete_indel = mainType[-1] + ':' + mainType_ins[0] + ':Ins:' + subType + ':' + mainType_ins[1]
														potential_sequence = ''
														for r in range(random_number, random_number+i,1):
															potential_sequence += tsb_ref[sequence[r]][1]
														seq_final = ''
														for r in range (random_number-1, i+random_number,1):
															seq_final += tsb_ref[sequence[r]][1]

														if vcf:
															print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
														else:	
															print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",tsb_ref[sequence[random_number-1]][1],".",seq_final,".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)

														if seqInfo:
															print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)

														if not overlap:
															for z in range (random_number-1, random_number+i+1,1):
																recorded_positions.add(z)

														# Updates the chromosome with the given mutation if desired
														if updating:
															sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Ins')
															location_range = len(sequence)

														# Remove one mutation count for the current INDEL and update all relevant data structures
														try:
															mutationsCountTSB[tsb_type][complete_indel] -= 1
															# print(chrom, tsb_type, mutationsCountTSB[tsb_type])
															u = 0
														except:
															sys.stderr.write(mainType_ins, complete_indel, indel_types[tsb_type], mutationsCountTSB)
															# print(mainType_ins, complete_indel, indel_types[tsb_type], mutationsCountTSB)
														if mutationsCountTSB[tsb_type][complete_indel] == 0:
															del mutationsCountTSB[tsb_type][complete_indel]
															indel_lengths[tsb_type].remove(int(mainType_ins[0]))
															repeat_lengths[tsb_type].remove(repeat_count+1+repeat_count_ins)
															if len(indel_types[tsb_type][mainType_ins]) > 1:
																del indel_types[tsb_type][mainType_ins][0]
															else:
																del indel_types[tsb_type][mainType_ins]
														break


												# Simulates microhomology deletions
												else:
													if repeat_count != 0 or i == 1:
														break
													max_repeat_length = max(repeat_lengths[tsb_type])+1
													if max_repeat_length > i:
														max_repeat_length = i


													# Counts homology in the forward direction
													homology_size1 = 0
													for k in range (1, 6, 1):
														seq = ''
														try:
															for r in range (random_number+i, random_number+k+i, 1):
																if r > len(sequence):
																	break
																seq += tsb_ref[sequence[r]][1]
														except:
															break
														if seq == inDel[:k]:
															homology_size1 += 1
														else:
															break

													# Counts homology in the reverse direction
													homology_size2  = 0
													for l in range (1, 6, 1):
														seq = ''
														try:
															for r in range (random_number-l, random_number, 1):
																if r > len(sequence):
																	break
																seq += tsb_ref[sequence[r]][1]
														except:
															break
														if seq == inDel[-l:]:
															homology_size2 += 1
														else:
															break

													# Assigns a new naming convetion for the current INDEL
													subType = 'M'
													if i > 5 and (homology_size1 >= 5 or homology_size2 >= 5):
														mainType1 = str(i-1) + str(homology_size1) + subType
														mainType2 = str(i-1) + str(homology_size2) + subType
													else:
														mainType1 = str(i) + str(homology_size1) + subType
														mainType2 = str(i) + str(homology_size2) + subType
													complete_indel = None
													
													mainType1 += bias
													mainType2 += bias

													# Checks to see if the forward homology is desired and that the reverse
													# homology equals 0. 
													if mainType1 in indel_types[tsb_type] and homology_size2 == 0:

															if not overlap:
																for p in range(random_number, random_number+i, 1):
																	if p in recorded_positions:
																		stop_flag = True
															if stop_flag:
																break

															# Reassigns the original naming convention for the INDEL and writes it to the output file
															# complete_indel = mainType[-1] + ':' + mainType1[0] + ':Del:' + subType + ':' + mainType1[1]
															complete_indel = bias + ':' + mainType1[0] + ':Del:' + subType + ':' + mainType1[1]

															seq_final = ''
															for r in range(random_number-1,i+random_number,1):
																seq_final += tsb_ref[sequence[r]][1]

															if vcf:
																print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
															else:
																print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",seq_final,".",tsb_ref[sequence[random_number-1]][1],".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
															
															if seqInfo:
																print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
															
															if not overlap:
																for z in range (random_number-1, random_number+i+1,1):
																	recorded_positions.add(z)

															# Updates the chromosome with the current INDEL if desired
															if updating:
																sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Del')
													 
															# Remove one mutation count for the current INDEL and update all relevant data structures       
															mutationsCountTSB[tsb_type][complete_indel] -= 1
															u = 0

															if mutationsCountTSB[tsb_type][complete_indel] == 0:
																del mutationsCountTSB[tsb_type][complete_indel]
																indel_lengths[tsb_type].remove(i)
																repeat_lengths[tsb_type].remove(homology_size1)
																if len(indel_types[tsb_type][mainType1]) > 1:
																	del indel_types[tsb_type][mainType1][0]
																else:
																	del indel_types[tsb_type][mainType1]
															break


													# Checks to see if the reverse homology is desired and that the forward
													# homology equals 0.
													elif mainType2 in indel_types[tsb_type] and homology_size1 == 0:
														seq = ''
														try:
															for r in range (random_number+i, random_number+(2*i),1):
																if r > len(sequence):
																	break
																seq += tsb_ref[sequence[r]][1]
														except:
															break

														seq2 = ''
														try:
															for r in range(random_number, random_number+i,1):
																if r > len(sequence):
																	break
																seq2 += tsb_ref[sequence[r]][1]
														except:
															break

														seq3 = ''
														try:
															for r in range(random_number-i, random_number, 1):
																if r > len(sequence):
																	break
																seq3 += tsb_ref[sequence[r]][1]
														except:
															break

														if indel_types[tsb_type][mainType2][0] == 'Del' and seq != seq2 and seq3 != seq2:
															if not overlap:
																for p in range(random_number, random_number+i, 1):
																	if p in recorded_positions:
																		stop_flag = True
															if stop_flag:
																break

															# Reassigns the original naming convention for the INDEL and writes it to the output file
															complete_indel = bias + ':' + mainType2[0] + ':Del:' + subType + ':' + mainType2[1]
															seq_final = ''
															for r in range (random_number-1,i+random_number,1):
																seq_final += tsb_ref[sequence[r]][1]

															if vcf:
																print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", seq_final,"\t",tsb_ref[sequence[random_number-1]][1],"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
															else:	
																print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",seq_final,".",tsb_ref[sequence[random_number-1]][1],".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
											
															if seqInfo:
																print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
															
															if not overlap:
																for z in range (random_number-1, random_number+i+1,1):
																	recorded_positions.add(z)

															# Updates the chromosome with the current INDEL if desired
															if updating:
																sequence = update_chromosome(sequence, random_number, sequence[random_number:i+random_number], 'Del')
																location_range = len(sequence)

															# Remove one mutation count for the current INDEL and update all relevant data structures
															mutationsCountTSB[tsb_type][complete_indel] -= 1
															u = 0
															if mutationsCountTSB[tsb_type][complete_indel] == 0:
																del mutationsCountTSB[tsb_type][complete_indel]
																indel_lengths[tsb_type].remove(i)
																repeat_lengths[tsb_type].remove(homology_size2)
																if len(indel_types[tsb_type][mainType2]) > 1:
																	del indel_types[tsb_type][mainType2][0]
																else:
																	del indel_types[tsb_type][mainType2]
															break




									# Simulates the insertion INDELs that have 0 repeats
									if not any(indel_types_O[tsb_type]):
										mutationsCountTSB[tsb_type] = {}
									while (any(indel_types_O[tsb_type]) == True):
										if u > 1000000:
											# logging.info(sample + " " + mutationsCount)
											# print (sample, mutationsCount, chrom)
											# print(indel_types_O[tsb_type])
											indel_types_O[tsb_type] = {}
											mutationsCountTSB[tsb_type] = {}
											u = 0
										# Randomly chooses a location on the current chromosome
										if not bed:
											if tsb_type != 'Q':
												location_range = chrom_bias_lengths[tsb_type][-1]
												random_range = fastrand.pcg32bounded(location_range)
												specific_range = bisect.bisect_left(chrom_bias_lengths[tsb_type], random_range)
												random_number = (chrom_bias_lengths[tsb_type][specific_range] - random_range) + chrom_bias[tsb_type][specific_range][0]
											else:
												random_number = fastrand.pcg32bounded(location_range)

										else:
											if tsb_type != 'Q':
												location_range = len(chrom_range)
												random_range = fastrand.pcg32bounded(location_range)
												random_number = chrom_range[random_range]
											else:
												random_number = fastrand.pcg32bounded(location_range)
												random_number = chrom_range[random_number]

										if random_number in recorded_positions and not overlap:
											continue										
										u += 1

										strand = "+"
										# Assigns a subtype for the chosen position on the chromosome
										if tsb_ref[sequence[random_number]][1] == 'T' or tsb_ref[sequence[random_number]][1] == 'A':
											subType = 'T'
											if tsb_ref[sequence[random_number]][1] == 'A':
												strand = "-"
										elif tsb_ref[sequence[random_number]][1] == 'G' or tsb_ref[sequence[random_number]][1] == 'C':
											subType = 'C'
											if tsb_ref[sequence[random_number]][1] == 'G':
												strand = "-"
										else:
											strand = "Q"
											break

										stop_flag = False
										for indels_O in indel_types_O[tsb_type]:

											# For INDELs of length 1, if the subType for the mutation does not match
											# the desired INDEL, break and find a new position
											if indels_O[0] == '1':
												if indels_O[2] != subType:
													break

											# # Assigns the original naming convention
											# complete_indel = indels_O[-1] + ':' + indels_O[0] + ':Ins:' + indels_O[2] + ':' + indels_O[1]

											bias = tsb_ref[sequence[random_number]][0]
											# Randomly assigns a base for insertions of length 1 based upon
											# the subType at the chosen position in the chromosome.
											potential_sequence = ''
											if int(indels_O[0]) == 1:
												if subType == 'T':
													potential_sequence += random_base(0,1)
												else:
													potential_sequence += random_base(2,3)

											# Randomly chooses bases until the insertion length equals the desired
											# INDEL length
											else:
												for m in range (0, int(indels_O[0]), 1):
													potential_sequence += random_base(0,3)

											if bias != 'N':
												if not bool(re.match("^[CT]*$", potential_sequence)) and not bool(re.match("^[GA]*$", potential_sequence)):
													bias = 'Q'

												else:
													if bool(re.match("^[GA]*$", potential_sequence)):
														bias = revbias(bias)	
											if tsb_type != bias:
												continue
											# Ensures that the bases are known and that there are no repeats around the insertion
											try:
												seq = ''
												for r in range (random_number-int(indels_O[0]),random_number,1):
													if r > len(sequence):
														break
													seq += tsb_ref[sequence[r]][1]
												seq2 = ''
												for r in range (random_number,random_number+int(indels_O[0]),1):
													if r > len(sequence):
														break
													seq2 += tsb_ref[sequence[r]][1]
											except:
												break

											if "N" not in seq:
												if seq2 != potential_sequence and tsb_ref[sequence[random_number-int(indels_O[0])]][1] != potential_sequence:
													if not overlap:
														for p in range(random_number, random_number+i, 1):
															if p in recorded_positions:
																stop_flag = True
													if stop_flag:
														break

													# Assigns the original naming convention
													complete_indel = bias + ':' + indels_O[0] + ':Ins:' + indels_O[2] + ':' + indels_O[1]

													if vcf:
														print (''.join([chrom,"\t",str(random_number),"\t",sample,"\t", tsb_ref[sequence[random_number-1]][1],"\t",tsb_ref[sequence[random_number-1]][1]+potential_sequence,"\t.\tSimulations\t",genome,"\t",complete_indel,"\t",strand]), file=out_vcf)
													else:	
														print ('\t'.join([".",".","sims",genome,chrom,str(random_number),str(random_number),strand,".","ID",tsb_ref[sequence[random_number-1]][1],".",tsb_ref[sequence[random_number-1]][1]+potential_sequence,".\t.",sample+"_"+str(simulations), complete_indel]), file=out_vcf)
													
													if seqInfo:
														print(''.join([sample, "\t",chrom,  "\t", str(random_number),  "\t",complete_indel, "\t", strand + "1"]), file=outSeq)
													if not overlap:
														for z in range (random_number-1, random_number+int(indels_O[0])+1,1):
															recorded_positions.add(z)

													# Remove one mutation count for the current INDEL and update all relevant data structures
													mutationsCountTSB[tsb_type][complete_indel] -= 1
													u = 0
													if mutationsCountTSB[tsb_type][complete_indel] == 0:
														del mutationsCountTSB[tsb_type][complete_indel]
														indel_lengths[tsb_type].remove(int(indels_O[0]))
														repeat_lengths[tsb_type].remove(int(indels_O[1]))
														del indel_types_O[tsb_type][indels_O]
													break

													# Updates the chromosome with the current INDEL if desired
													if updating:
														sequence = update_chromosome(sequence, random_number, sequence[random_number:int(indels_O[0])+random_number], 'Ins')
														location_range = len(sequence)
											break

			
									
									
						# Simulates 96, 1536, and DINUCs        
						elif sim == 2 or sim == 3 or sim == 5 or sim == 1:

							# Organizes nucleotide keys for later reference.
							for nuc in nuc_keys:
								if sim == 1:
									base_keys.append(nuc[0])
								elif sim == 2:
									base_keys.append(nuc[0] + nuc[2] + nuc[6])
								elif sim == 3:
									base_keys.append(nuc[0:2] + nuc[3] + nuc[7:])
								elif sim == 5:
									base_keys.append(nuc[0:2])
							
							# Simulates until all mutations have been assigned.
							l = 0            
							while (any(mutationsCount) == True):

								# Picks a random location to throw a mutation limited to the
								random_number = fastrand.pcg32bounded(location_range)
								if bed:
									random_number = chrom_range[random_number]
								

								# If a specific mutation cannot be assinged after x iterations,
								# skip that nucleotide context. Helps to prevent the simulation from
								# stalling on a rare/non-existent mutation
								l += 1
								if l > 1000000:
									logging.info(sample + " " + mutationsCount)
									print (sample, mutationsCount, chrom)
								# 	if sample not in left_over_mutations:
								# 		left_over_mutations[sample] = {}
								# 		left_over_mutations[sample][simulations] = {context:None}
								# 		left_over_mutations[sample][simulations][context] = mutationsCount
								# 	else:
								# 		#if simulations not in left_over_mutations[sample].keys():
								# 		if simulations not in left_over_mutations[sample]:
								# 			left_over_mutations[sample][simulations] = {context:mutationsCount}
								# 		else:
								# 			left_over_mutations[sample][simulations][context] = mutationsCount
									mutationsCount = {}

								# For DINUC context: organizes nucleotide references
								if sim == 5:
									mutNuc = ''.join([tsb_ref[base][1] for base in sequence[random_number:random_number+2]])
									revCompMutNuc = revcompl(mutNuc)       

								# For all other contexts: organize nucleotide references                 
								else:
									mutNuc = ''.join([tsb_ref[base][1] for base in sequence[random_number - mut_start:random_number + mut_start+1]])
									revCompMutNuc = revcompl(mutNuc)
					
								
								# If the nucleotide is desired (present in the mutation dictionary), write
								# it to the output file and update the dictionary
								bases = None
								if mutNuc in base_keys:
									nucIndex = base_keys.index(mutNuc)
									if nuc_keys[nucIndex] in mutationsCount and mutationsCount[nuc_keys[nucIndex]] != 0:
										
										# Exclude mutations if new position overlaps an existing mutation or if it is within
										# the user-specified spacing (default=1bp to exclude DBSs)
										if not overlap:
											if random_number in recorded_positions:
												continue
											mnv_flag = False
											for i in range(random_number - spacing, random_number + spacing, 1):
												if i in recorded_positions:
													mnv_flag = True
													break
											if mnv_flag:
												continue



										if sim != 5 and sim != 4 and sim != 7:
											context_up = 'SNP'
											bases = nuc_keys[nucIndex][mut_save+2]

											if vcf:
												print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", nuc_keys[nucIndex][mut_save],"\t",nuc_keys[nucIndex][mut_save+2],"\t.\tSimulations\t",genome,"\t",mutNuc,"\t","+1"]), file=out_vcf)
											else:
												print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"+1",".","SNP", nuc_keys[nucIndex][mut_save],".",nuc_keys[nucIndex][mut_save+2],".\t.",sample+"_"+str(simulations), mutNuc]), file=out_vcf)
											if seqInfo:
												mutNuc_seq = ''.join([tsb_ref[base][1] for base in sequence[random_number - 2:random_number + 3]])
												print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",mutNuc_seq[0:2],"[",nuc_keys[nucIndex][mut_save],">",nuc_keys[nucIndex][mut_save+2],"]",mutNuc_seq[3:], "\t","+1"]), file=outSeq)
						
											recorded_positions.add(random_number)



										elif sim == 5:
											context_up = 'DBS'
											bases = nuc_keys[nucIndex][mut_save+3:mut_save+5]
											if vcf:
												print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", nuc_keys[nucIndex][mut_save:mut_save+2],"\t",nuc_keys[nucIndex][mut_save+3:mut_save+5],"\t.\tSimulations\t",genome,"\t",mutNuc,"\t","+1"]), file=out_vcf)
											else:	
												print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"+1",".","DBS", nuc_keys[nucIndex][mut_save:mut_save+2],".",nuc_keys[nucIndex][mut_save+3:mut_save+5],".\t.",sample+"_"+str(simulations), mutNuc]), file=out_vcf)
											if seqInfo:
												print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",nuc_keys[nucIndex][mut_save:mut_save+2],">",nuc_keys[nucIndex][mut_save+3:mut_save+5],"\t","+1"]), file=outSeq)										
											recorded_positions.add(random_number)

										mutationsCount[nuc_keys[nucIndex]] -= 1
										l = 0
										if mutationsCount[nuc_keys[nucIndex]] == 0:
											del mutationsCount[nuc_keys[nucIndex]]
											del nuc_keys[nucIndex]
											del base_keys[nucIndex] 
						
								# If the reverse complement of the nucleotide context is desired,
								# write it to the output file as the reverse complement.
								elif revCompMutNuc in base_keys:
									nucIndex = base_keys.index(revCompMutNuc)
									if nuc_keys[nucIndex] in mutationsCount and mutationsCount[nuc_keys[nucIndex]] != 0:

										# Exclude mutations if new position overlaps an existing mutation or if it is within
										# the user-specified spacing (default=1bp to exclude DBSs)
										if not overlap:
											if random_number in recorded_positions:
												continue
											mnv_flag = False
											for i in range(random_number - spacing, random_number + spacing, 1):
												if i in recorded_positions:
													mnv_flag = True
													break
											if mnv_flag:
												continue


										if sim != 5 and sim != 4 and sim != 7:
											context_up = 'SNP'
											bases = revcompl(nuc_keys[nucIndex][mut_save+2])
											if vcf:
												print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", revcompl(nuc_keys[nucIndex][mut_save]),"\t",revcompl(nuc_keys[nucIndex][mut_save+2]),"\t.\tSimulations\t",genome,"\t",revcompl(revCompMutNuc),"\t","-1"]), file=out_vcf)
											else:
												print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"-1",".","SNP", revcompl(nuc_keys[nucIndex][mut_save]),".",revcompl(nuc_keys[nucIndex][mut_save+2]),".\t.",sample+"_"+str(simulations), revcompl(revCompMutNuc)]), file=out_vcf)
											if seqInfo:
												revCompMutNuc_seq = revcompl(''.join([tsb_ref[base][1] for base in sequence[random_number - 2:random_number + 3]]))
												print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",revCompMutNuc_seq[0:2],"[",nuc_keys[nucIndex][mut_save],">",nuc_keys[nucIndex][mut_save+2],"]",revCompMutNuc_seq[3:], "\t","-1"]), file=outSeq)
											recorded_positions.add(random_number)
										elif sim == 5:
											context_up = 'DBS'
											bases = revcompl(nuc_keys[nucIndex][mut_save+3:mut_save+5])
											if vcf:
												print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", revcompl(nuc_keys[nucIndex][mut_save:mut_save+2]),"\t",revcompl(nuc_keys[nucIndex][mut_save+3:mut_save+5]),"\t.\tSimulations\t",genome,"\t",revcompl(revCompMutNuc),"\t","-1"]), file=out_vcf)
											else:
												print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"-1",".","DBS", revcompl(nuc_keys[nucIndex][mut_save:mut_save+2]),".",revcompl(nuc_keys[nucIndex][mut_save+3:mut_save+5]),".\t.",sample+"_"+str(simulations), revcompl(revCompMutNuc)]), file=out_vcf)
							
											if seqInfo:
												print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",nuc_keys[nucIndex][mut_save:mut_save+2],">",nuc_keys[nucIndex][mut_save+3:mut_save+5], "\t","-1"]), file=outSeq)
											recorded_positions.add(random_number)
										mutationsCount[nuc_keys[nucIndex]] -= 1
										l = 0 
										if mutationsCount[nuc_keys[nucIndex]] == 0:
											del mutationsCount[nuc_keys[nucIndex]]
											del nuc_keys[nucIndex]
											del base_keys[nucIndex]

								# If the user specified udpating mutations, proceeds to update the chromosome
								# with the current mutation
								# if updating and bases != None:
								# 	bias = tsb_ref[sequence[random_number]][0]
								# 	new_bases = ''.join([str(tsb_ref_rev[bias][base]) for base in bases])
								# 	sequence = update_chromosome(sequence, random_number, new_bases, context_up)
								# 	# location_range = len(sequence)
								if updating and bases != None:
									bias = tsb_ref[sequence[random_number]][0]
									new_bases = ''.join([str(tsb_ref_rev[bias][base]) for base in bases])
									for i in range (0, len(new_bases), 1):
										sequence[random_number+i] = int(new_bases[i])
									location_range = len(sequence)

						# Simulates TSB [24, 384, 6144]
						else:
							tsb = ['TU','B','N']
							mutationsCountTSB = {'TU':{},'B':{},'N':{}}
							base_keys = {'T':[],'U':[],'B':[],'N':[]}
							nuc_keys_tsb = {'T':[],'U':[],'B':[],'N':[]}
							if sim == 8:
								mutationsCountTSB['Q'] = {}
								tsb.append('Q')
								base_keys['Q'] = []
								nuc_keys_tsb['Q'] = []



							# Organizes nucleotide keys for later reference.
							for nuc in nuc_keys:
								if nuc[0] == 'T' or nuc[0] == 'U':
									mutationsCountTSB['TU'][nuc] = mutationsCount[nuc]
								else:
									mutationsCountTSB[nuc[0]][nuc] = mutationsCount[nuc]
								bias = nuc[0]
								nuc_keys_tsb[bias].append(nuc)
								if sim == 10:
									base_keys[nuc[0]].append(nuc[0] + nuc[2])
								elif sim == 4:
									base_keys[nuc[0]].append(nuc[0] + nuc[2] + nuc[4] + nuc[8])
								elif sim == 7:
									base_keys[nuc[0]].append(nuc[0] + nuc[2:4] + nuc[5] + nuc[9:])
								elif sim == 8:
									base_keys[nuc[0]].append(nuc[0] + nuc[2:4])

							# Simulates until all mutations have been assigned.
							for tsb_type in tsb:
								l = 0 
								while (any(mutationsCountTSB[tsb_type]) == True):

									# Picks a random location to throw a mutation limited to the
									# length of the current chromosome
									if not bed:
										if tsb_type != 'Q':
											location_range = chrom_bias_lengths[tsb_type][-1]
											random_range = fastrand.pcg32bounded(location_range)
											specific_range = bisect.bisect_left(chrom_bias_lengths[tsb_type], random_range)
											random_number = (chrom_bias_lengths[tsb_type][specific_range] - random_range) + chrom_bias[tsb_type][specific_range][0]
										else:
											location_range = len(sequence)
											random_number = fastrand.pcg32bounded(location_range)

									else:
										if tsb_type != 'Q':
											location_range = len(chrom_range)
											random_range = fastrand.pcg32bounded(location_range)
											random_number = chrom_range[random_range]
										else:
											random_number = fastrand.pcg32bounded(location_range)
											random_number = chrom_range[random_number]
										

									# If a specific mutation cannot be assinged after x iterations,
									# skip that nucleotide context. Helps to prevent the simulation from
									# stalling on a rare/non-existent mutation
									l += 1
									if l > 1000000:
										print(mutationsCountTSB[tsb_type])
										if sample not in left_over_mutations:
											left_over_mutations[sample] = {}
											left_over_mutations[sample][simulations] = {context:{}}
											for nuc in mutationsCountTSB[tsb_type]:
												left_over_mutations[sample][simulations][context][nuc] = mutationsCountTSB[tsb_type][nuc]
										else:
											if simulations not in left_over_mutations[sample]:
												left_over_mutations[sample][simulations] = {}
												left_over_mutations[sample][simulations] = {context:{}}
												for nuc in mutationsCountTSB[tsb_type]:
													left_over_mutations[sample][simulations][context][nuc] = mutationsCountTSB[tsb_type][nuc]
											else:
												for nuc in mutationsCountTSB[tsb_type]:
													left_over_mutations[sample][simulations][context][nuc] = mutationsCountTSB[tsb_type][nuc]

										mutationsCountTSB[tsb_type] = {}
										l = 0
						
									# Only for TSB simulations: organizes nucleotide references
									# nuc_bias = tsb_type
									nuc_bias = tsb_ref[sequence[random_number]][0]

									if sim == 8:
										mutNuc = ''.join([tsb_ref[base][1] for base in sequence[random_number:random_number+2]])
									else:
										mutNuc = ''.join([tsb_ref[base][1] for base in sequence[random_number - mut_start:random_number + mut_start+1]])
									
									mutNuc = nuc_bias + mutNuc
									revCompMutNuc = revbias(nuc_bias) + revcompl(mutNuc[1:])  

									if tsb_type == 'Q':
										if mutNuc in dinuc_non_tsb or revCompMutNuc in dinuc_non_tsb:
											l -= 1
											continue	

									# If the nucleotide is desired (present in the mutation dictionary), write
									# it to the output file and update the dictionary
									bases = None
									if mutNuc in base_keys[nuc_bias]:
										if random_range not in ranges:
											ranges[random_range] = 1
										else:
											ranges[random_range] += 1
										nucIndex = base_keys[nuc_bias].index(mutNuc)
										if nuc_keys_tsb[nuc_bias][nucIndex] in mutationsCountTSB[tsb_type] and mutationsCountTSB[tsb_type][nuc_keys_tsb[nuc_bias][nucIndex]] != 0:		

											# Exclude mutations if new position overlaps an existing mutation or if it is within
											# the user-specified spacing (default=1bp to exclude DBSs)
											if not overlap:
												if random_number in recorded_positions:
													continue
												mnv_flag = False
												for i in range(random_number - spacing, random_number + spacing, 1):
													if i in recorded_positions:
														mnv_flag = True
														break
												if mnv_flag:
													continue	

											
											if sim == 8:
												bases = nuc_keys_tsb[nuc_bias][nucIndex][mut_save+3:mut_save+5]
												context_up = 'DBS'
												if vcf:
													print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", nuc_keys_tsb[nuc_bias][nucIndex][mut_save:mut_save+2],"\t",nuc_keys_tsb[nuc_bias][nucIndex][mut_save+3:mut_save+5],"\t.\tSimulations\t",genome,"\t",mutNuc,"\t","+1"]), file=out_vcf)
												
												else:
													print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"+1",".","DBS", nuc_keys_tsb[nuc_bias][nucIndex][2:4],".",nuc_keys_tsb[nuc_bias][nucIndex][5:],".\t.",sample+"_"+str(simulations), mutNuc]), file=out_vcf)
												if seqInfo:
													print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",nuc_keys_tsb[nuc_bias][nucIndex][2:4],">",nuc_keys_tsb[nuc_bias][nucIndex][5:],"\t","+1"]), file=outSeq)										
												recorded_positions.add(random_number)
											else:
												bases = nuc_keys_tsb[nuc_bias][nucIndex][mut_save+2]
												context_up = 'SNP'
												if vcf:
													print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", nuc_keys_tsb[nuc_bias][nucIndex][mut_save],"\t",nuc_keys_tsb[nuc_bias][nucIndex][mut_save+2],"\t.\tSimulations\t",genome,"\t",mutNuc,"\t","+1"]), file=out_vcf)
												else:
													print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"+1",".","SNP", nuc_keys_tsb[nuc_bias][nucIndex][mut_save],".",nuc_keys_tsb[nuc_bias][nucIndex][mut_save+2],".\t.",sample+"_"+str(simulations), mutNuc]), file=out_vcf)

												if seqInfo:
													mutNuc_seq = ''.join([tsb_ref[base][1] for base in sequence[random_number - 2:random_number + 3]])
													mutNuc_seq = nuc_bias + ":" + mutNuc_seq
													print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",mutNuc_seq[0:4], "[",nuc_keys_tsb[nuc_bias][nucIndex][mut_save],">",nuc_keys_tsb[nuc_bias][nucIndex][mut_save+2],"]",mutNuc_seq[5:], "\t","+1"]), file=outSeq)
										
											recorded_positions.add(random_number)
											mutationsCountTSB[tsb_type][nuc_keys_tsb[nuc_bias][nucIndex]] -= 1

											l = 0
							
											if mutationsCountTSB[tsb_type][nuc_keys_tsb[nuc_bias][nucIndex]] == 0:
												del mutationsCountTSB[tsb_type][nuc_keys_tsb[nuc_bias][nucIndex]]
												del nuc_keys_tsb[nuc_bias][nucIndex]
												del base_keys[nuc_bias][nucIndex] 

									# If the reverse complement of the nucleotide context is desired,
									# write it to the output file as the reverse complement.
									elif revCompMutNuc in base_keys[revbias(nuc_bias)]:
										if random_range not in ranges:
											ranges[random_range] = 1
										else:
											ranges[random_range] += 1
										rev_tsb_type = revbias(nuc_bias)

										nucIndex = base_keys[rev_tsb_type].index(revCompMutNuc)
										if nuc_keys_tsb[rev_tsb_type][nucIndex] in mutationsCountTSB[tsb_type] and mutationsCountTSB[tsb_type][nuc_keys_tsb[rev_tsb_type][nucIndex]] != 0:
										
											# Exclude mutations if new position overlaps an existing mutation or if it is within
											# the user-specified spacing (default=1bp to exclude DBSs)
											if not overlap:
												if random_number in recorded_positions:
													continue
												mnv_flag = False
												for i in range(random_number - spacing, random_number + spacing, 1):
													if i in recorded_positions:
														mnv_flag = True
														break
												if mnv_flag:
													continue	

											if sim == 8:
												context_up = 'DBS'
												bases = revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save+3:mut_save+5])

												if vcf:
													print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save:mut_save+2]),"\t",revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save+3:mut_save+5]),"\t.\tSimulations\t",genome,"\t",revCompMutNuc,"\t","-1"]), file=out_vcf)

												else:
													print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"-1",".","DBS", revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][2:4]),".",revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][5:]),".\t.",sample+"_"+str(simulations), revCompMutNuc]), file=out_vcf)
												if seqInfo:
													print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",nuc_keys_tsb[rev_tsb_type][nucIndex][2:4],">",nuc_keys_tsb[rev_tsb_type][nucIndex][5:], "\t","-1"]), file=outSeq)
											
											else:
												context_up = 'SNP'
												bases = revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save+2])
												if vcf:
													print (''.join([chrom,"\t",str(random_number+1),"\t",sample,"\t", revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save]),"\t",revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save+2]),"\t.\tSimulations\t",genome,"\t",revCompMutNuc,"\t","-1"]), file=out_vcf)
												else:
													print ('\t'.join([".",".","sims",genome,chrom,str(random_number+1),str(random_number+1),"-1",".","SNP", revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save]),".",revcompl(nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save+2]),".\t.",sample+"_"+str(simulations), revCompMutNuc]), file=out_vcf)

												if seqInfo:
													revCompMutNuc_seq = revbias(nuc_bias) + ":" + revcompl(''.join([tsb_ref[base][1] for base in sequence[random_number - 2:random_number + 3]]))
													print(''.join([sample, "\t",chrom,  "\t", str(random_number+1),  "\t",revCompMutNuc_seq[0:4],"[",nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save],">",nuc_keys_tsb[rev_tsb_type][nucIndex][mut_save+2],"]",revCompMutNuc_seq[5:], "\t","-1"]), file=outSeq)
											
											
											recorded_positions.add(random_number)
											mutationsCountTSB[tsb_type][nuc_keys_tsb[rev_tsb_type][nucIndex]] -= 1
											l = 0

											if mutationsCountTSB[tsb_type][nuc_keys_tsb[rev_tsb_type][nucIndex]] == 0:
												del mutationsCountTSB[tsb_type][nuc_keys_tsb[rev_tsb_type][nucIndex]]
												del nuc_keys_tsb[rev_tsb_type][nucIndex]
												del base_keys[rev_tsb_type][nucIndex]

									# If the user specified udpating mutations, proceeds to update the chromosome
									# with the current mutation
									if updating and bases != None:
										bias = tsb_ref[sequence[random_number]][0]
										new_bases = ''.join([tsb_ref_rev[bias][base] for base in bases])
										sequence = update_chromosome(sequence, random_number, new_bases, context_up)
										location_range = len(sequence)
							

				simulations -= 1
				if seqInfo:
					outSeq.flush()
					outSeq.close()

		# with open("/Users/ebergstr/Desktop/" + chrom + "_tsbStats.txt", "w") as out_stats:
		# 	print(len(ranges), file=out_stats)
		# 	for x, y in ranges.items():
		# 		if y > 1:
		# 			print(str(x) + "\t" + str(y) + "\t" + str(chrom_range[x]), file=out_stats)

			# a = [x + ":" + str(y) for x,y in ranges.items() if y >1]
			# print(a, chrom, len(a), )

		out_log = open(log_file, 'a')
		print("Chromosome " + chrom + " done", file=out_log)
		out_log.close()
		print("         Chromosome " + chrom + " done")

	return(chrom)
	