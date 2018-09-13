#!/usr/bin/env python3
import pickle
import os

chromosome_path = '/Users/ebergstr/Desktop/Perl_tests/testCode/simulation_code_python/mutation_simulation/references/chromosomes/tsb/mm10/'
chromosome_BED_path = '/Users/ebergstr/Desktop/Perl_tests/testCode/simulation_code_python/mutation_simulation/references/chromosomes/tsb_BED/mm10/'
chromosomes = ['X', 'Y', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
			   '13', '14', '15', '16', '17', '18', '19']#, '20', '21', '22']


if not os.path.exists(chromosome_BED_path):
	os.makedirs(chromosome_BED_path)

for chrom in chromosomes:
	with open (chromosome_path+chrom+"_192.txt", "rb") as f, open(chromosome_BED_path + chrom + "_BED_TSB.txt", 'w') as out:
		print("<CHROM>\t<START>\t<END>\t<TSB>", file=out)
		chrom_tsb = f.read()
		first_tsb = chrom_tsb[0]
		current_range = [0]
		for i in range (1, len(chrom_tsb), 1):
			if chrom_tsb[i] != first_tsb:
				current_range.append(i-1)
				print(chrom + "\t" + str(current_range[0]) + "\t" + str(current_range[1]) + "\t" + str(chrom_tsb[i-1]), file=out)
				first_tsb = chrom_tsb[i]
				current_range = [i]
			else:
				continue
	print("chromosome ", chrom, "done")
