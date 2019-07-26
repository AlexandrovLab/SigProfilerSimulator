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
import multiprocessing as mp
import numpy as np
import platform
import SigProfilerSimulator as sigSim
import SigProfilerMatrixGenerator as sig
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGenerator as matRef
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from . import mutational_simulator as simScript
from SigProfilerMatrixGenerator.scripts import save_context_distribution as context_dist



def SigProfilerSimulator (project, project_path, genome, contexts, exome=None, simulations=1, updating=False, bed_file=None, overlap=False, gender='female', seqInfo=False, chrom_based=False, seed_file=None):
	'''
	contexts -> [] must be a list
	'''
	start_run = time.time()

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
	error_file = project_path + 'logs/SigProfilerSimulator_' + project + "_" + genome + "_" + str(time_stamp) + ".err"
	log_file = project_path + 'logs/SigProfilerSimulator_' + project + "_" + genome + "_" + str(time_stamp) + ".out"

	if not os.path.exists(project_path + "logs/"):
		os.makedirs(project_path + "logs/")

	if os.path.exists(error_file):
		os.system("rm " + error_file)
	if os.path.exists(log_file):
		os.system("rm " + log_file)


	sys.stderr = open(error_file, 'w')
	log_out = open(log_file, 'w')
	log_out.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
	log_out.write("-------System Info-------\n")
	log_out.write("Operating System Name: "+ platform.uname()[0]+"\n"+"Nodename: "+ platform.uname()[1]+"\n"+"Release: "+ platform.uname()[2]+"\n"+"Version: "+ platform.uname()[3]+"\n")
	log_out.write("\n-------Python and Package Versions------- \n")
	log_out.write("Python Version: "+str(platform.sys.version_info.major)+"."+str(platform.sys.version_info.minor)+"."+str(platform.sys.version_info.micro)+"\n")
	log_out.write("SigProfilerSimulator Version: "+sigSim.__version__+"\n")
	log_out.write("SigProfilerMatrixGenerator Version: "+sig.__version__+"\n")
	log_out.write("numpy version: "+np.__version__+"\n")
	
	log_out.write("\n-------Vital Parameters Used for the execution -------\n")
	log_out.write("Project: {}\nGenome: {}\nInput File Path: {}\ncontexts: {}\nexome: {}\nsimulations: {}\nupdating: {}\nbed_file: {}\noverlap: {}\ngender: {}\nseqInfo: {}\nchrom_based: {}\nseed_file: {}\n".format(project, project_path, genome, contexts, str(exome), str(simulations),  str(updating), str(bed_file), str(overlap), gender, str(seqInfo), str(chrom_based), str(seed_file)))
	log_out.write("\n-------Date and Time Data------- \n")
	tic = datetime.datetime.now()
	log_out.write("Date and Clock time when the execution started: "+str(tic)+"\n\n\n")
	


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
	for context in contexts:
		matrix_path = project_path + "output/"
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
		if chrom_based:
			if os.path.exists (catalogue_file + '.chr1') == False:
				if os.path.exists (vcf_files_2) == False and len(os.listdir(vcf_files_1)) == 0:
					print ("Please place your vcf files for each sample into the 'references/vcf_files/[project]/' directory. Once you have done that, rerun this script.")
				else:
					print("Matrices per chromosomes do not exist. Creating the matrix files now.")
					matGen.SigProfilerMatrixGeneratorFunc(project, genome, project_path ,plot=False, exome=exome, bed_file=bed_file, chrom_based=True)
					print("The matrix file has been created. Continuing with simulations...")
			if os.path.exists (catalogue_file) == False:
				if os.path.exists (vcf_files_2) == False and len(os.listdir(vcf_files_1)) == 0:
					print ("Please place your vcf files for each sample into the 'references/vcf_files/[project]/' directory. Once you have done that, rerun this script.")
				else:
					print(catalogue_file + " does not exist. Creating the matrix file now.")
					matGen.SigProfilerMatrixGeneratorFunc(project, genome, project_path ,plot=False, exome=exome, bed_file=bed_file)
					print("The matrix file has been created. Continuing with simulations...")


		else:	
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
	elif exome:
		output_path = project_path + "output/simulations/" + project + '_simulations_' + genome + '_' + context_string + '_exome/'
	else:
		output_path = project_path + "output/simulations/" + project + '_simulations_' + genome + '_' + context_string + '/'

	if os.path.exists(output_path):
		shutil.rmtree(output_path)
		os.makedirs(output_path)

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
	if chrom_based:
		sample_names, mut_prep, mut_dict = simScript.mutation_preparation_chromosomes(catalogue_files, matrix_path, chromosomes, project, log_file)
		reference_sample = sample_names[0]
	else:
		sample_names, mut_prep = simScript.mutation_preparation(catalogue_files, log_file)
		reference_sample = sample_names[0]
		mut_dict = simScript.mut_tracker(sample_names,  mut_prep, reference_sample, nucleotide_context_files, chromosome_string_path, genome, chromosomes, bed_file, log_file)
	
	# Set-up parallelization:
	processors = mp.cpu_count()
	max_seed = processors
	if processors > len(chromosomes):
		max_seed = len(chromosomes)
	pool = mp.Pool(max_seed)

	chrom_break = len(chromosomes)/max_seed
	chromosomes_parallel = [[] for i in range(max_seed)]

	chrom_bin = 0
	for chrom in chromosomes:
		if chrom_bin == max_seed:
			chrom_bin = 0
		chromosomes_parallel[chrom_bin].append(chrom)
		chrom_bin += 1

	# Generate unique seeds for each process
	log_out.write("\n-------Seeds for random number generation per process------- \n")
	seeds = []
	if seed_file == None:
		ref_dir, tail = os.path.split(os.path.dirname(os.path.abspath(__file__)))
		seed_file = ref_dir + "/SigProfilerSimulator/seeds.txt"
	with open(seed_file) as f:
		for i in range (0, max_seed, 1):
			new_seed = int(f.readline().strip())
			seeds.append(new_seed)
			log_out.write("Process " + str(i) + ": " + str(new_seed) + "\n")

	log_out.write("\n\n\n-------Runtime Checkpoints------- \n")
	log_out.close()

	# For chroms in chromosomes_parallel:
	for i in range (0, len(chromosomes_parallel), 1):
		pool.apply_async(simScript.simulator, args=(sample_names,  mut_prep, mut_dict, chromosome_string_path, tsb_ref, tsb_ref_rev, simulations, seeds[i], output_path, updating, chromosomes_parallel[i], project, genome, bed, bed_file, contexts, exome, overlap, project_path, seqInfo, log_file))
	pool.close()
	pool.join()



	for sample in sample_names:
		concat_files = os.listdir(output_path + sample + "/")
		for i in range (1, simulations + 1, 1):
			with open(output_path + sample + "/" + sample + "_" + str(i) + ".vcf", "wb") as f:
				for chrom in chromosomes:
					with open(output_path + sample + "/" + sample + "_" + str(i) + "_" + chrom + ".vcf",'rb') as fd:
						shutil.copyfileobj(fd, f)
					os.remove(output_path + sample + "/" + sample + "_" + str(i) + "_" + chrom + ".vcf")


	end_run = time.time()
	run_time = end_run - start_run
	log_out = open(log_file, 'a')
	print("Simulation completed\nJob took " , run_time, " seconds", file=log_out)
	print("Simulation completed\nJob took " , run_time, " seconds")
	log_out.close()





