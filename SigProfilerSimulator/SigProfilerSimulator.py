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
import logging
import datetime
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as matRef
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
#import mutational_simulator_with_range_option_simulate_multiple_at_once_sig_option_192_alternate_io_buffering_new_update_bi_Y as simScript
from . import mutational_simulator as simScript
from SigProfilerMatrixGenerator.scripts import save_context_distribution as context_dist
#import save_context_distribution_96_192_1536_DINUC_with_range_option_bi_fix as context_dist

start_run = time.time()


def SigProfilerSimulator (project, project_path, genome, contexts, exome=None, simulations=1, updating=False, bed_file=None, Signatures=False, overlap=False, gender='male', seqInfo=False):
	'''
	contexts -> [] must be a list
	'''

	# Ensures proper string for the project's path
	if project_path[-1] != "/":
		project_path += "/"

	# Sorts the user-provided contexts
	contexts.sort(reverse=True)

	bed = False
	if bed_file:
		bed = True
	exome_file = exome

	# Asigns a species based on the genome parameter
	species = None
	if genome.upper() == 'GRCH37' or genome.upper() == 'GRCH38': 
		species = "homo_sapiens"
	elif genome.upper() == 'MM10' or genome.upper() == 'MM9': 
		species = "mus_musculus"
	else:
		print(genome + " is not supported. The following genomes are supported:\nGRCh37, GRCh38, mm9, mm10.")


	############################## References ###########################################################################################################
	chromosomes = ['X', 'Y', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
				   '13', '14', '15', '16', '17', '18', '19', '20', '21', '22']
	
	tsb_ref = {0:['N','A'], 1:['N','C'], 2:['N','G'], 3:['N','T'],
			   4:['T','A'], 5:['T','C'], 6:['T','G'], 7:['T','T'],
			   8:['U','A'], 9:['U','C'], 10:['U','G'], 11:['U','T'],
			   12:['B','A'], 13:['B','C'], 14:['B','G'], 15:['B','T'],
			   16:['N','N'], 17:['T','N'], 18:['U','N'], 19:['B','N']}

	tsb_ref_rev = {'N':{'A':0, 'C':1, 'G':2, 'T':3, 'N':16},
			   	   'T':{'A':4, 'C':5, 'G':6, 'T':7, 'N':17},
			       'U':{'A':8, 'C':9, 'G':10, 'T':11, 'N':18},
			       'B':{'A':12, 'C':13, 'G':14, 'T':15, 'N':19}}
			       
	if species == 'mus_musculus':
		chromosomes = chromosomes[:21]

	if gender == 'female' or gender.upper() == 'F':
		chromosomes.remove('Y')

	
	
	############################## Log and Error Files ##################################################################################################
	time_stamp = datetime.date.today()
	error_file = project_path + 'logs/sigProfilerSimulator_' + project + "_" + genome + "_" + str(time_stamp) + ".err"
	log_file = project_path + 'logs/sigProfilerSimulator_' + project + "_" + genome + "_" + str(time_stamp) + ".out"

	if not os.path.exists(project_path + "logs/"):
		os.makedirs(project_path + "logs/")

	if os.path.exists(error_file):
		os.system("rm " + error_file)
	if os.path.exists(log_file):
		os.system("rm " + log_file)

	sys.stderr = open(error_file, 'w')
	logging.basicConfig(filename=log_file, level=logging.INFO)



	############################## Pre-simulation Checks ##################################################################################################
	# Ensures that the chromosome strings are saves properly:
	chromosome_string_path, ref_dir = matRef.reference_paths(genome)
	if os.path.exists(chromosome_string_path) == False or len(os.listdir(chromosome_string_path)) <= len(chromosomes):
		print("The chromosome strings were not saved properly or have not been created yet. Rerun the SigProfilerMatrixGenerator isntall script.")

	# Ensures that the chromosome proportions are saved: 
	if os.path.exists(chromosome_string_path + genome + "_proportions.txt") == False:
		print("Chromosome proportion file does not exist. Creating now...")
		chromosomeProbs = simScript.chrom_proportions(chromosome_string_path, genome, chromosomes)
		print("Chromosome proportion file created. Proceeding with simulation...")

	if bed_file:
		print("Creating a chromosome proportion file for the given BED file ranges...")
		chromosomeProbs = simScript.chrom_proportions_BED(bed_file, chromosome_string_path, genome, chromosomes)

	# Ensures that the mutational matrices exist:
	catalogue_files = {}
	
	#matrix_path = "references/matrix/" + project + "/"
	for context in contexts:
		matrix_path = project_path + "output/"
		#matrix_path = matrix_path + context + "/"
		if context == 'DINUC' or context == 'DBS':
			context_folder = 'DBS'
			matrix_path = matrix_path + context_folder + "/"
			file_name = ".DBS78"
		elif context == 'INDEL' or context == 'ID':
			context_folder = 'ID'
			matrix_path = matrix_path + context_folder + "/"
			file_name = '.ID83'
		else:
			context_folder = 'SBS'
			matrix_path = matrix_path + context_folder + "/"
			file_name = '.SBS' + context
		if exome:
			catalogue_file = matrix_path + project + file_name + '.exome'
		else:
			if bed_file:
				catalogue_file = matrix_path + project + file_name + '.region'
			else:
				catalogue_file = matrix_path + project + file_name + '.all'
	
		catalogue_files[context] = catalogue_file

		vcf_files_1 = project_path
		vcf_files_2 = project_path + "input/"
		parent_dir = os.getcwd()
		matrix_dir = "scripts/"
		if os.path.exists (catalogue_file) == False:
			if os.path.exists (vcf_files_2) == False and len(os.listdir(vcf_files_1)) == 0:
				print ("Please place your vcf files for each sample into the 'references/vcf_files/[project]/' directory. Once you have done that, rerun this script.")
			else:
				print(catalogue_file + " does not exist. Creating the matrix file now.")
				matGen.SigProfilerMatrixGeneratorFunc(project, genome, project_path ,plot=False, exome=exome, bed_file=bed_file)
				print("The matrix file has been created. Continuing with simulations...")


	# Esnures that the nucleotide context files are saved properly
	nucleotide_context_files = {}
	for context in contexts:
		nucleotide_context_file = chromosome_string_path.split("/")
		ref_path = nucleotide_context_file[:-3]
		ref_path = '/'.join([x for x in ref_path])
		nucleotide_context_file = ref_path + '/context_distributions/'
		
		if bed_file:
			nucleotide_context_file += "context_distribution_" + genome + "_" + context + "_" + gender + "_BED.csv"
		else:
			if exome:
				nucleotide_context_file += "context_distribution_" + genome + "_" + context + "_" + gender + "_exome.csv"
			else:
				nucleotide_context_file += "context_distribution_" + genome + "_" + context + "_" + gender + ".csv"

		nucleotide_context_files[context] = nucleotide_context_file

		if os.path.exists(nucleotide_context_file) == False and (context != 'INDEL' and context != 'ID'):
			print("The context distribution file does not exist. This file needs to be created before simulating. This may take several hours...")
			if bed:
				output_file = ref_path + 'context_distributions/context_distribution_' + genome + "_" + context + "_" + gender + '_BED.csv'
				context_dist.context_distribution_BED(context, output_file, chromosome_string_path, chromosomes, bed, bed_file, exome, exome_file, genome, ref_path, tsb_ref)
			elif exome:
				output_file = ref_path + 'context_distributions/context_distribution_' + genome + "_" + context + "_" + gender + '_exome.csv'
				context_dist.context_distribution_BED(context, output_file, chromosome_string_path, chromosomes, bed, bed_file, exome, exome_file, genome, ref_dir, tsb_ref)
			else:
				output_file = ref_path + 'context_distributions/context_distribution_' + genome + "_" + context + "_" + gender + '.csv'
				context_dist.context_distribution(context, output_file, chromosome_string_path, chromosomes, tsb_ref)
			print("Context distribution file successfully created. Proceeding with simulation...")



	############################## Set-up output files ##################################################################################################
	context_string = "_".join(contexts)
	if bed_file:
		output_path = project_path + "output/simulations/" + project + '_simulations_' + genome + '_' + context_string + '_BED/'
	elif Signatures:
		output_path = project_path + "output/simulations/" + project + '_simulations_' + genome + '_' + context_string + '_signature_based/'
	elif exome:
		output_path = project_path + "output/simulations/" + project + '_simulations_' + genome + '_' + context_string + '_exome/'
	else:
		output_path = project_path + "output/simulations/" + project + '_simulations_' + genome + '_' + context_string + '/'



	############################## Set parameters for simulation ##################################################################################################

	sim = None
	mut_start = None
	mut_length = None
	if context == '96':
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
	elif context == '3072':
		sim = 7
		mut_save = 5
		mut_start = 2



	############################## Begin the simulation process ##################################################################################################
	if Signatures:
		mut_prep = simScript.mutation_preparation_sig(catalogue_files)
	else:
		mut_prep = simScript.mutation_preparation(catalogue_files)
	reference_sample = mut_prep[0][0]
	mut_dict = simScript.mut_tracker(mut_prep[0], mut_prep[1], reference_sample, nucleotide_context_files, chromosome_string_path, genome, chromosomes, bed_file)
	simScript.simulator(mut_prep[0], mut_prep[1], mut_dict, chromosome_string_path, tsb_ref, tsb_ref_rev, simulations, output_path, updating, chromosomes, project, genome, bed, bed_file, contexts, exome, overlap, project_path, seqInfo)
	end_run = time.time()
	run_time = end_run - start_run
	logging.info("Simulation completed\nJob took " + str(run_time) + " seconds")
	print("Simulation completed\nJob took " , run_time, " seconds")




