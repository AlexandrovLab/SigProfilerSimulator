#!/usr/bin/env python3

# Author: Erik Bergstrom

# Contact: ebergstr@eng.ucsd.edu


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
import pandas as pd
import SigProfilerSimulator as sigSim
import SigProfilerMatrixGenerator as sig
from SigProfilerMatrixGenerator.scripts import MutationMatrixGenerator as matRef
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
from . import mutational_simulator as simScript
from SigProfilerMatrixGenerator.scripts import save_context_distribution as context_dist


def context_identifier(mutation):
    SBS_prefix = set([0, 1, 2, 3, 4])
    SBS_ref = set([1, 3, 4, 5, 6])
    DBS_prefix = set([0, 1, 3, 5, 7])
    DBS_ref = set([2, 4, 6])
    mutation_prefix = len(mutation.split("[")[0])
    mutation_ref = len(mutation.split(">")[0])
    tsb_length = len(mutation.split(":")[0])
    mutation_length = len(mutation)
    context = ""
    nuc = mutation
    nuc = nuc.replace("[", "")
    nuc = nuc.replace("]", "")

    if mutation_prefix in SBS_prefix and mutation_ref in SBS_ref:
        # context = "SBS"
        context = ""
        nuc = nuc.split(">")[0] + nuc.split(">")[1][1:]
        if mutation_length == 3:
            context += "6"
        elif mutation_length == 5:
            context += "24"
        elif mutation_length == 7:
            context += "96"
        elif mutation_length == 9:
            if mutation_ref == 5:
                context += "384"
            else:
                context += "1536"
        elif mutation_length == 11:
            context += "6144"
        else:
            print(mutation, " is not supported by this function.")
            sys.exit()
    elif mutation_prefix in DBS_prefix and mutation_ref in DBS_ref:
        context = "DBS"
        nuc = nuc.split(">")[0]
        if mutation_length == 5:
            context += 78
        elif mutation_length == 7:
            context += "186"
        else:
            print(mutation, " is not supported by this function")
    else:
        print(mutation, " is not supported by this function.")
        sys.exit()

    return (context, nuc)


def probability(
    chromosome=None,
    position=None,
    mutation=None,
    context=None,
    genome=None,
    mutation_count=1,
    mutation_file=None,
    exome=False,
):

    chromosome_string_path, ref_dir = matRef.reference_paths(genome)
    if not mutation_file:
        if not genome:
            print("No genome provided")
            sys.exit()
        if not chromosome:
            print("No chromosome provided")
            sys.exit()
        if not position:
            print("No position provided")
            sys.exit()
        if not mutation:
            print("No mutation provided.")
            sys.exit()

        context, nuc = context_identifier(mutation)
        if exome:
            context += "_exome"
        nucleotide_context_file = (
            ref_dir
            + "/references/chromosomes/context_distributions/"
            + "context_counts_"
            + genome
            + "_"
            + context
            + ".csv"
        )
        count_mat = pd.read_csv(
            nucleotide_context_file, sep=",", header=0, index_col=[0]
        )

        nucleotide_count = count_mat.loc[nuc, chromosome]

        prob = mutation_count / nucleotide_count
        print(
            "The probabilty of seeing",
            mutation,
            " on chromosome",
            chromosome,
            "at position",
            str(position),
            "is equal to:\n\n\t\t",
            str(prob),
        )

        # else:
        # 	first_line = True
        # 	output_path = os.path.dirname(mutation_file)
        # 	with open(mutation_file) as f, open(output_path + "probabilties.txt", "w") as out:
        # 		for lines in f:
        # 			lines = lines.strip().split()
        # 			chrom = lines[0]
        # 			pos = lines[1]
        # 			mutation = lines[2]
        # 			exome = lines[3]
        # 			context = context_identifier(mutation)
        # 			if first_line:
        # 				nucleotide_context_file += "context_counts_" + genome + "_" + context + ".csv"
        # 				count_mat = pd.DataFrame.from_csv(nucleotide_context_file, sep=',', header=0)
        # 				first_line = False
        # 			prob = mutation_count/nucleotide_count
        # 			print("\t".join([chrom, pos, mutation, prob]), file=out)

        pass


def SigProfilerSimulator(
    project,
    project_path,
    genome,
    contexts,
    exome=None,
    simulations=1,
    updating=False,
    bed_file=None,
    overlap=False,
    gender="female",
    seqInfo=False,
    chrom_based=False,
    seed_file=None,
    spacing=1,
    noisePoisson=False,
    noiseUniform=0,
    cushion=100,
    region=None,
    vcf=False,
    mask=None,
):
    """
    contexts -> [] must be a list
    """
    print(
        "\n======================================\n        SigProfilerSimulator        \n======================================\n\nChecking for all reference files and relevant matrices..."
    )
    start_run = time.time()

    # Ensures proper string for the project's path
    if project_path[-1] != "/":
        project_path += "/"

    # Sorts the user-provided contexts
    contexts.sort(reverse=True)

    bed = False
    if bed_file:
        bed = True
    exome_file = None

    # Asigns a species based on the genome parameter
    species = None
    if (
        genome.upper() == "GRCH37"
        or genome.upper() == "GRCH38"
        or "GRCH37" in genome.upper()
        or "GRCH38" in genome.upper()
    ):
        species = "homo_sapiens"
    elif (
        genome.upper() == "MM10"
        or genome.upper() == "MM9"
        or "MM9" in genome.upper()
        or "MM10" in genome.upper()
        or genome.upper() == "MM39"
        or "MM39" in genome.upper()
    ):
        species = "mus_musculus"
    else:
        species = "custom"

    ############################## References ###########################################################################################################
    chromosomes = [
        "X",
        "Y",
        "1",
        "2",
        "3",
        "4",
        "5",
        "6",
        "7",
        "8",
        "9",
        "10",
        "11",
        "12",
        "13",
        "14",
        "15",
        "16",
        "17",
        "18",
        "19",
        "20",
        "21",
        "22",
    ]

    tsb_ref = {
        0: ["N", "A"],
        1: ["N", "C"],
        2: ["N", "G"],
        3: ["N", "T"],
        4: ["T", "A"],
        5: ["T", "C"],
        6: ["T", "G"],
        7: ["T", "T"],
        8: ["U", "A"],
        9: ["U", "C"],
        10: ["U", "G"],
        11: ["U", "T"],
        12: ["B", "A"],
        13: ["B", "C"],
        14: ["B", "G"],
        15: ["B", "T"],
        16: ["N", "N"],
        17: ["T", "N"],
        18: ["U", "N"],
        19: ["B", "N"],
    }

    tsb_ref_rev = {
        "N": {"A": 0, "C": 1, "G": 2, "T": 3, "N": 16},
        "T": {"A": 4, "C": 5, "G": 6, "T": 7, "N": 17},
        "U": {"A": 8, "C": 9, "G": 10, "T": 11, "N": 18},
        "B": {"A": 12, "C": 13, "G": 14, "T": 15, "N": 19},
    }

    if species == "mus_musculus":
        chromosomes = chromosomes[:21]

    chromosome_string_path, ref_dir = matRef.reference_paths(genome)
    if species == "custom":
        chromosome_string_path, ref_dir = matRef.reference_paths(genome)
        chromosomes = os.listdir(chromosome_string_path)
        if ".DS_Store" in chromosomes:
            chromosomes.remove(".DS_Store")

        chromosomes = [x.split(".")[0] for x in chromosomes if len(x.split(".")[0]) < 8]
        if genome == "yeast":
            chromosomes = sorted(
                chromosomes,
                key=lambda x: (
                    [
                        "I",
                        "II",
                        "III",
                        "IV",
                        "V",
                        "VI",
                        "VII",
                        "VIII",
                        "IX",
                        "X",
                        "XI",
                        "XII",
                        "XIII",
                        "XIV",
                        "XV",
                        "XVI",
                    ].index(x)
                ),
            )
    if gender == "female" or gender.upper() == "FEMALE":
        if "Y" in chromosomes:
            chromosomes.remove("Y")

    if region:
        chromosomes = [region]
    ############################## Log and Error Files ##################################################################################################
    time_stamp = datetime.date.today()
    error_file = (
        project_path
        + "logs/SigProfilerSimulator_"
        + project
        + "_"
        + genome
        + "_"
        + str(time_stamp)
        + ".err"
    )
    log_file = (
        project_path
        + "logs/SigProfilerSimulator_"
        + project
        + "_"
        + genome
        + "_"
        + str(time_stamp)
        + ".out"
    )

    if not os.path.exists(project_path + "logs/"):
        os.makedirs(project_path + "logs/")

    if os.path.exists(error_file):
        # os.system("rm " + error_file)
        os.remove(error_file)
    if os.path.exists(log_file):
        # os.system("rm " + log_file)
        os.remove(log_file)

    sys.stderr = open(error_file, "w")
    log_out = open(log_file, "w")
    log_out.write("THIS FILE CONTAINS THE METADATA ABOUT SYSTEM AND RUNTIME\n\n\n")
    log_out.write("-------System Info-------\n")
    log_out.write(
        "Operating System Name: "
        + platform.uname()[0]
        + "\n"
        + "Nodename: "
        + platform.uname()[1]
        + "\n"
        + "Release: "
        + platform.uname()[2]
        + "\n"
        + "Version: "
        + platform.uname()[3]
        + "\n"
    )
    log_out.write("\n-------Python and Package Versions------- \n")
    log_out.write(
        "Python Version: "
        + str(platform.sys.version_info.major)
        + "."
        + str(platform.sys.version_info.minor)
        + "."
        + str(platform.sys.version_info.micro)
        + "\n"
    )
    log_out.write("SigProfilerSimulator Version: " + sigSim.__version__ + "\n")
    log_out.write("SigProfilerMatrixGenerator Version: " + sig.__version__ + "\n")
    log_out.write("numpy version: " + np.__version__ + "\n")

    log_out.write("\n-------Vital Parameters Used for the execution -------\n")
    log_out.write(
        "Project: {}\nGenome: {}\nInput File Path: {}\ncontexts: {}\nexome: {}\nsimulations: {}\nupdating: {}\nbed_file: {}\noverlap: {}\ngender: {}\nseqInfo: {}\nchrom_based: {}\nseed_file: {}\n".format(
            project,
            genome,
            project_path,
            contexts,
            str(exome),
            str(simulations),
            str(updating),
            str(bed_file),
            str(overlap),
            gender,
            str(seqInfo),
            str(chrom_based),
            str(seed_file),
        )
    )
    log_out.write("\n-------Date and Time Data------- \n")
    tic = datetime.datetime.now()
    log_out.write(
        "Date and Clock time when the execution started: " + str(tic) + "\n\n\n"
    )

    ############################## Pre-simulation Checks ##################################################################################################
    # Ensures that the chromosome strings are saves properly:
    chromosome_string_path, ref_dir = matRef.reference_paths(genome)
    if os.path.exists(chromosome_string_path) == False or len(
        os.listdir(chromosome_string_path)
    ) < len(chromosomes):
        print(
            "     The chromosome strings were not saved properly or have not been created yet. Please refer to the SigProfilerMatrixGenerator README for installation instructions:\n\thttps://github.com/AlexandrovLab/SigProfilerMatrixGenerator"
        )
        sys.exit()
    # Ensures that the chromosome proportions are saved:
    if os.path.exists(chromosome_string_path + genome + "_proportions.txt") == False:
        print("     Chromosome proportion file does not exist. Creating now...", end="")
        chromosomeProbs = simScript.chrom_proportions(
            chromosome_string_path, genome, chromosomes
        )
        print("Completed!")

    if bed_file:
        print(
            "     Creating a chromosome proportion file for the given BED file ranges...",
            end="",
        )
        chromosomeProbs = simScript.chrom_proportions_BED(
            bed_file, chromosome_string_path, genome, chromosomes
        )
        print("Completed!")

    # Ensures that the mutational matrices exist:
    catalogue_files = {}
    for context in contexts:
        matrix_path = project_path + "output/"
        if context == "DINUC" or "DBS" in context:
            context_folder = "DBS"
            matrix_path = matrix_path + context_folder + "/"
            if context == "DBS" or context == "DINUC" or context == "78":
                file_name = ".DBS78"
            else:
                file_name = "." + context
        elif context == "INDEL" or "ID" in context or "415" in context:
            context_folder = "ID"
            matrix_path = matrix_path + context_folder + "/"
            if context == "INDEL" or context == "ID" or context == "83":
                file_name = ".ID83"
            else:
                file_name = "." + context
        else:
            context_folder = "SBS"
            matrix_path = matrix_path + context_folder + "/"
            file_name = ".SBS" + context

        if exome:
            catalogue_file = matrix_path + project + file_name + ".exome"
        else:
            if bed_file:
                catalogue_file = matrix_path + project + file_name + ".region"
            else:
                catalogue_file = matrix_path + project + file_name + ".all"

        catalogue_files[context] = catalogue_file

        vcf_files_1 = project_path
        vcf_files_2 = project_path + "input/"
        parent_dir = os.getcwd()
        matrix_dir = "scripts/"
        if chrom_based:
            if os.path.exists(catalogue_file + ".chr1") == False:
                if (
                    os.path.exists(vcf_files_2) == False
                    and len(os.listdir(vcf_files_1)) == 0
                ):
                    print(
                        "     Please place your vcf files for each sample into the 'references/vcf_files/[project]/' directory. Once you have done that, rerun this script."
                    )
                else:
                    print(
                        "     Matrices per chromosomes do not exist. Creating the matrix files now."
                    )
                    matGen.SigProfilerMatrixGeneratorFunc(
                        project,
                        genome,
                        project_path,
                        plot=False,
                        exome=exome,
                        bed_file=bed_file,
                        chrom_based=True,
                        cushion=cushion,
                    )
                    # print("The matrix file has been created. Continuing with simulations...")
            if os.path.exists(catalogue_file) == False:
                if (
                    os.path.exists(vcf_files_2) == False
                    and len(os.listdir(vcf_files_1)) == 0
                ):
                    print(
                        "     Please place your vcf files for each sample into the 'references/vcf_files/[project]/' directory. Once you have done that, rerun this script."
                    )
                else:
                    print(
                        "     "
                        + catalogue_file
                        + " does not exist. Creating the matrix file now."
                    )
                    matGen.SigProfilerMatrixGeneratorFunc(
                        project,
                        genome,
                        project_path,
                        plot=False,
                        exome=exome,
                        bed_file=bed_file,
                        cushion=cushion,
                    )
                    # print("The matrix file has been created. Continuing with simulations...")

        else:
            if os.path.exists(catalogue_file) == False:  # or bed_file:
                if (
                    os.path.exists(vcf_files_2) == False
                    and len(os.listdir(vcf_files_1)) == 0
                ):
                    print(
                        "     Please place your vcf files for each sample into the 'references/vcf_files/[project]/' directory. Once you have done that, rerun this script."
                    )
                else:
                    print(
                        "     "
                        + catalogue_file
                        + " does not exist. Creating the matrix file now."
                    )
                    matGen.SigProfilerMatrixGeneratorFunc(
                        project,
                        genome,
                        project_path,
                        plot=False,
                        exome=exome,
                        bed_file=bed_file,
                        cushion=cushion,
                    )
                    # print("The matrix file has been created. Continuing with simulations...")

    if exome:
        exome_file = (
            ref_dir
            + "/references/chromosomes/exome/"
            + genome
            + "/"
            + genome
            + "_exome.interval_list"
        )

    # Esnures that the nucleotide context files are saved properly
    nucleotide_context_files = {}
    for context in contexts:
        nucleotide_context_file = chromosome_string_path.split("/")
        ref_path = nucleotide_context_file[:-3]
        ref_path = "/".join([x for x in ref_path])
        nucleotide_context_file = ref_path + "/context_distributions/"

        # genome_original = genome
        # if 'havana' in genome:
        # 	genome = genome.split("_")[0]

        if bed_file:
            if region:
                nucleotide_context_file += (
                    "context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + ".csv"
                )
            else:
                nucleotide_context_file += (
                    "context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + "_BED.csv"
                )
        else:
            if exome:
                nucleotide_context_file += (
                    "context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + "_exome.csv"
                )
            else:
                nucleotide_context_file += (
                    "context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + ".csv"
                )

        if context == "288":
            nucleotide_context_file = nucleotide_context_file.split("_")
            nucleotide_context_file[5] = "384"
            nucleotide_context_file = "_".join([x for x in nucleotide_context_file])
        elif context == "4608":
            nucleotide_context_file = nucleotide_context_file.split("_")
            nucleotide_context_file[4] = "6144"
            nucleotide_context_file = "_".join([x for x in nucleotide_context_file])
        nucleotide_context_files[context] = nucleotide_context_file
        if os.path.exists(nucleotide_context_file) == True and bed and not region:
            os.remove(nucleotide_context_file)

        if os.path.exists(nucleotide_context_file) == False and (
            context != "INDEL" and context != "ID" and context != "ID415"
        ):
            print(
                "     The context distribution file does not exist. This file needs to be created before simulating. This may take several hours..."
            )
            if bed:
                output_file = (
                    ref_path
                    + "/context_distributions/context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + "_BED.csv"
                )
                context_dist.context_distribution_BED(
                    context,
                    output_file,
                    chromosome_string_path,
                    chromosomes,
                    bed,
                    bed_file,
                    exome,
                    exome_file,
                    genome,
                    ref_path,
                    tsb_ref,
                    gender,
                )
            elif exome:
                output_file = (
                    ref_path
                    + "/context_distributions/context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + "_exome.csv"
                )
                context_dist.context_distribution_BED(
                    context,
                    output_file,
                    chromosome_string_path,
                    chromosomes,
                    bed,
                    bed_file,
                    exome,
                    exome_file,
                    genome,
                    ref_dir,
                    tsb_ref,
                    gender,
                )
            else:
                output_file = (
                    ref_path
                    + "/context_distributions/context_distribution_"
                    + genome
                    + "_"
                    + context
                    + "_"
                    + gender
                    + ".csv"
                )
                context_dist.context_distribution(
                    context,
                    output_file,
                    chromosome_string_path,
                    chromosomes,
                    tsb_ref,
                    genome,
                )
            print("     The context distribution file has been created!")
            if gender == "female" or gender.upper() == "FEMALE":
                if "Y" in chromosomes:
                    chromosomes.remove("Y")

    ############################## Set-up output files ##################################################################################################
    context_string = "_".join(contexts)
    if bed_file:
        output_path = (
            project_path
            + "output/simulations/"
            + project
            + "_simulations_"
            + genome
            + "_"
            + context_string
            + "_BED/"
        )
    elif exome:
        output_path = (
            project_path
            + "output/simulations/"
            + project
            + "_simulations_"
            + genome
            + "_"
            + context_string
            + "_exome/"
        )
    else:
        output_path = (
            project_path
            + "output/simulations/"
            + project
            + "_simulations_"
            + genome
            + "_"
            + context_string
            + "/"
        )

    if os.path.exists(output_path):
        shutil.rmtree(output_path)
        os.makedirs(output_path)
    else:
        os.makedirs(output_path)

    if "M" in chromosomes:
        chromosomes.remove("M")
    if "MT" in chromosomes:
        chromosomes.remove("MT")
    ############################## Begin the simulation process ##################################################################################################
    print()
    if chrom_based:
        sample_names, mut_prep, mut_dict = simScript.mutation_preparation_chromosomes(
            catalogue_files, matrix_path, chromosomes, project, log_file
        )
        reference_sample = sample_names[0]
    elif region:
        sample_names, mut_prep, mut_dict = simScript.mutation_preparation_region(
            catalogue_files, matrix_path, project, log_file, region
        )
        reference_sample = sample_names[0]
    else:
        sample_names, mut_prep = simScript.mutation_preparation(
            catalogue_files, log_file
        )
        reference_sample = sample_names[0]
        mut_dict = simScript.mut_tracker(
            sample_names,
            mut_prep,
            reference_sample,
            nucleotide_context_files,
            chromosome_string_path,
            genome,
            chromosomes,
            bed_file,
            log_file,
        )

    if vcf:
        if "" in sample_names:
            sample_names.remove("")
        for sample in sample_names:
            if not os.path.exists(output_path + sample + "/"):
                os.makedirs(output_path + sample + "/")

    # Set-up parallelization:
    processors = mp.cpu_count()
    max_seed = processors
    if processors > len(chromosomes):
        max_seed = len(chromosomes)
    pool = mp.Pool(max_seed)

    chrom_break = len(chromosomes) / max_seed
    chromosomes_parallel = [[] for i in range(max_seed)]

    chrom_bin = 0
    for chrom in chromosomes:
        if chrom_bin == max_seed:
            chrom_bin = 0
        chromosomes_parallel[chrom_bin].append(chrom)
        chrom_bin += 1

    iterations_parallel = [[] for i in range(max_seed)]
    iter_bin = 0
    for i in range(1, simulations + 1, 1):
        if iter_bin == max_seed:
            iter_bin = 0
        iterations_parallel[iter_bin].append(i)
        iter_bin += 1

    # Generate unique seeds for each process
    log_out.write("\n-------Seeds for random number generation per process------- \n")
    seeds = []
    if seed_file == None:
        ref_dir, tail = os.path.split(os.path.dirname(os.path.abspath(__file__)))
        seed_file = ref_dir + "/SigProfilerSimulator/seeds.txt"
    with open(seed_file) as f:
        for i in range(0, max_seed, 1):
            new_seed = int(int(f.readline().strip()) / time.time())
            seeds.append(new_seed)
            log_out.write("Process " + str(i) + ": " + str(new_seed) + "\n")

    log_out.write("\n\n\n-------Runtime Checkpoints------- \n")
    log_out.close()

    if exome:
        bed = True
        bed_file = (
            ref_dir
            + "/SigProfilerMatrixGenerator/references/chromosomes/exome/"
            + genome
            + "/"
            + genome
            + "_exome.interval_list"
        )

    if seqInfo:
        seqOut_path = project_path + "output/vcf_files/simulations/"
        if not os.path.exists(seqOut_path):
            os.makedirs(seqOut_path)

        for context in contexts:
            if not os.path.exists(seqOut_path + context + "/"):
                os.makedirs(seqOut_path + context + "/")
            else:
                print(seqOut_path + context + "/")
                shutil.rmtree(seqOut_path + context + "/")
                os.makedirs(seqOut_path + context + "/")

    pool = mp.Pool(max_seed)
    results = []
    for i in range(0, len(chromosomes_parallel), 1):
        mut_dict_parallel = {
            k1: {
                k2: {
                    k3: {
                        k4: v4 for k4, v4 in v3.items() if k4 in chromosomes_parallel[i]
                    }
                    for k3, v3 in v2.items()
                }
                for k2, v2 in v1.items()
            }
            for k1, v1 in mut_dict.items()
        }
        r = pool.apply_async(
            simScript.simulator,
            args=(
                sample_names,
                mut_dict_parallel,
                chromosome_string_path,
                tsb_ref,
                tsb_ref_rev,
                simulations,
                seeds[i],
                cushion,
                output_path,
                updating,
                chromosomes_parallel[i],
                project,
                genome,
                bed,
                bed_file,
                contexts,
                overlap,
                project_path,
                seqInfo,
                log_file,
                spacing,
                noisePoisson,
                noiseUniform,
                vcf,
                mask,
            ),
        )
        results.append(r)
    pool.close()
    pool.join()
    # simScript.simulator(sample_names, mut_dict, chromosome_string_path, tsb_ref, tsb_ref_rev, simulations, seeds[0], output_path, updating, chromosomes, project, genome, bed, bed_file, contexts, overlap, project_path, seqInfo, log_file, spacing, noisePoisson, noiseAWGN)
    for r in results:
        r.wait()
        if not r.successful():
            # Raises an error when not successful
            r.get()

    pool = mp.Pool(max_seed)

    # if region:
    bed = False

    for i in range(0, len(iterations_parallel), 1):
        r = pool.apply_async(
            simScript.combine_simulation_files,
            args=(
                iterations_parallel[i],
                output_path,
                chromosomes,
                sample_names,
                bed,
                exome,
                vcf,
            ),
        )
    pool.close()
    pool.join()

    for r in results:
        r.wait()
        if not r.successful():
            # Raises an error when not successful
            r.get()

    end_run = time.time()
    run_time = end_run - start_run
    log_out = open(log_file, "a")
    print("Simulation completed\nJob took ", run_time, " seconds", file=log_out)
    print("Simulation completed\nJob took ", run_time, " seconds")
    log_out.close()
    sys.stderr.close()
