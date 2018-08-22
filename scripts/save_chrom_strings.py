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

import os
import re


def save_chrom_strings ():
	current_dir = os.getcwd()
	ref_dir = re.sub('\/scripts$', '', current_dir)
	chrom_fasta_path = ref_dir + '/references/chromosomes/fasta/' + genome + "/"
	chrom_string_path = ref_dir + '/references/chromosomes/chrom_string/' + genome + '/'
	chrom_fasta_files = os.listdir(chrom_fasta_path)

	for files in chrom_fasta_files:
		file_name = files.split(".")
		chromosome = file_name[-2]
		if file_name[-4] == 'dna':
			with open(files) as chrom, open(chrom_string_path + chromosome + ".txt", 'w') as out:
				chromosome_final = ''
				for lines in chrom:
					line = lines.strip()
					chromosome_final += line


				print(chromosome_final, file=out)
				print("The string file for Chromosome " + chromosome + " has been created.")