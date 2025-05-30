[![Docs](https://img.shields.io/badge/docs-latest-blue.svg)](https://osf.io/usxjz/wiki/home/) [![License](https://img.shields.io/badge/License-BSD\%202--Clause-orange.svg)](https://opensource.org/licenses/BSD-2-Clause) [![Build Status](https://app.travis-ci.com/AlexandrovLab/SigProfilerSimulator.svg?branch=master)](https://app.travis-ci.com/AlexandrovLab/SigProfilerSimulator)

# SigProfilerSimulator 
SigProfilerSimulator allows realistic simulations of mutational signatures in cancer genomes. The tool can be used to simulate signatures of single point mutations, double point mutations, and insertion/deletions. Further, the tool makes use of SigProfilerMatrixGenerator and SigProfilerPlotting.  

# INTRODUCTION
The purpose of this document is to provide a guide for using the SigProfilerSimulator for simulating mutational signatures in cancer. This tool allows for realistic simulations of single point mutations, double point mutations, and insertions/deletions with the goal of providing a background model for statistical analysis. The simulations are performed in an unbiased fashion, relying on random chance as the main distribution and can be performed across the entire genome or limited to user-provided ranges. This tool currently supports the GRCh37, GRCh38, mm9, and mm10 assemblies, however, additional genomes may be installed. In addition, this tool makes use of SigProfilerMatrixGenerator and SigProfilerPlotting. An extensive Wiki page detailing the usage of this tool can be found at https://osf.io/usxjz/wiki/home/.

For users that prefer working in an R environment, a wrapper package is provided and can be found and installed from: https://github.com/AlexandrovLab/SigProfilerSimulatorR

![schematic](Figure1.png)

# PREREQUISITES
The framework is written in PYTHON, however, it also requires the following additional software with the given versions (or newer) and access to BASH:

-PYTHON (version 3.4 or newer)

-FASTRAND (Python module: https://github.com/lemire/fastrand/blob/master/README.md )

-SigProfilerMatrixGenerator (Current version: https://github.com/AlexandrovLab/SigProfilerMatrixGenerator )

-Desired reference genome (Follow the installation process on the SigProfilerMatrixGenerator README )

While the code was developed for application on a local computer, access to a cluster with greater computational power may be required for simulating a large number of mutations/samples.

# QUICK START GUIDE
This section will guide you through the minimum steps required to begin simulating mutations:
1. First, install the python package using pip. The R wrapper still requires the python package:
```
                          pip install SigProfilerSimulator
```
2. Place your vcf files in your desired output folder. It is recommended that you name this folder based on your project's name

3. From within a Python3 session, you can now simulate mutational patterns/signatures as follows:
```
$ python3
>> from SigProfilerSimulator import SigProfilerSimulator as sigSim
>> sigSim.SigProfilerSimulator("BRCA", "/Users/ebergstr/Desktop/BRCA/", "GRCh37", contexts=["96"], exome=None, simulations=100, updating=False, bed_file=None, overlap=False, gender='female',  chrom_based=False, seed_file=None, noisePoisson=False, noiseAWGN=0, cushion=100, region=None, vcf=False)
```
  The layout of the required parameters are as follows:
  
      SigProfilerSimulator ( project, project_path, genome, contexts)
            
  where project, project_path, and genome must be strings (surrounded by quotation marks, ex: "test"), and contexts is a list of the desired contexts to simulate (ex: contexts=["96", "ID"]) Optional
  parameters include:
  
### Parameters

| Category         | Parameter       | Type     | Description |
|------------------|----------------|----------|-------------|
| **Required**     | `project`      | `str`    | Name for the simulation run. Used as a prefix for output files. |
|                  | `project_path` | `str`    | Directory where output will be saved. Must contain an `input/` subdirectory with VCF or MAF files. |
|                  | `genome`       | `str`    | Reference genome to use (`"GRCh37"`, `"GRCh38"`, `"mm10"`, etc.). |
|                  | `contexts`     | `list`   | Mutation contexts to simulate (e.g., `["96"]`, `["ID"]`, `["96", "ID"]`). |
| **Optional**     | `exome`        | `bool`   | Simulates mutations on the exome of the reference genome. |
|                  | `simulations`  | `int`    | Number of iterations to simulate per sample. Default is `1`. |
|                  | `updating`     | `bool`   | Updates the chromosome with each mutation. Default is `False`. |
|                  | `bed_file`     | `str`    | Full path to a BED file for simulating on custom genomic regions. |
|                  | `overlap`      | `bool`   | Allows overlapping mutations along the chromosome. Default is `False`. |
|                  | `gender`       | `str`    | Gender of the simulated genome (`"male"` or `"female"`). Default is `"female"`. |
|                  | `chrom_based`  | `bool`   | maintains the same catalogs of mutations on a per chromosome basis. Default is `False`. |
|                  | `seed_file`    | `str`    | Path to a file containing a single master seed. This seed is used to generate a sequence of internal seeds for parallel simulation. If `None`, a new seed is generated and saved to `<project_path>/output/Simulator_seeds.txt`. ***To reproduce runs exactly, run with the same number of CPU cores.*** |
|                  | `noisePoisson` | `bool`   | Adds Poisson noise to the simulations. Default is `False`. |
|                  | `noiseAWGN`    | `float`  | Adds additive white Gaussian noise (±%) to the mutation counts per category (e.g., `noiseAWGN=5` allows ±2.5% variation). Default is `0`. |
|                  | `cushion`      | `int`    | Base-pair padding around exons or target regions during simulation. Default is `100`. |
|                  | `region`       | `str`    | Path to a target region panel for simulating on custom-defined regions. Default is whole-genome simulations. |
|                  | `vcf`          | `bool`   | If `True`, outputs simulated samples as individual VCF files per iteration per sample. Default is `False` (MAF format). |
|                  | `mask`         | `str`    | Path to a probability mask file (TSV format) with required columns: `Chromosome`, `Start`, `End`, `Probability`. Not compatible with `bed_file=True`. |


**INPUT FILE FORMAT**

This tool currently supports maf, vcf, simple text file, and ICGC formats. The user must provide variant data adhering to one of these four formats. If the users' files are in vcf format, each sample must be saved as a separate files. For an example input file, please download the simple text file "example.txt" from the following link: [example.txt](https://osf.io/7rdkh/)


**Output File Structure**

The output structure is divided into three folders: input, output, and logs. The input folder contains copies of the user-provided input files. The output folder contains
a DBS, SBS, ID, and simulations folder. The matrices are saved into the appropriate folders, and the simulations are found within a project specific folder under simulations. The logs folder contains the error and log files for the submitted job.


**SUPPORTED GENOMES**

This tool currently supports the following genomes:

GRCh38.p12 [GRCh38] (Genome Reference Consortium Human Reference 37), INSDC
Assembly GCA_000001405.27, Dec 2013. Released July 2014. Last updated January 2018. This genome was downloaded from ENSEMBL database version 93.38.

GRCh37.p13 [GRCh37] (Genome Reference Consortium Human Reference 37), INSDC
Assembly GCA_000001405.14, Feb 2009. Released April 2011. Last updated September 2013. This genome was downloaded from ENSEMBL database version 93.37. 

GRCm38.p6 [mm10] (Genome Reference Consortium Mouse Reference 38), INDSDC
Assembly GCA_000001635.8, Jan 2012. Released July 2012. Last updated March 2018. This genome was downloaded from ENSEMBL database version 93.38. 

GRCm37 [mm9] (Release 67, NCBIM37), INDSDC Assembly GCA_000001635.18.
Released Jan 2011. Last updated March 2012. This genome was downloaded from ENSEMBL database version release 67.

rn6 (Rnor_6.0) INSDC Assembly GCA_000001895.4, Jul 2014. Released Jun 2015. Last updated Jan 2017. 
This genome was downloaded from ENSEMBL database version 96.6.

yeast (Saccharomyces cerevisiae S288C; assembly R64-2-1). Released Nov 2014. 

**LOG FILES**

All errors and progress checkpoints are saved into *sigProfilerSimulator_[project]_[genome].err* and *sigProfilerSimulator_[project]_[genome].out*, respectively. 
For all errors, please email the error and progress log files to the primary contact under CONTACT INFORMATION.

**CITATION**

Erik N. Bergstrom, Mark Barnes, I�igo Martincorena, Ludmil B. Alexandrov
bioRxiv 2020.02.13.948422; doi: https://doi.org/10.1101/2020.02.13.948422
https://www.biorxiv.org/content/10.1101/2020.02.13.948422v1

**COPYRIGHT**

Copyright (c) 2020, Erik Bergstrom [Alexandrov Lab] All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 

**CONTACT INFORMATION**

Please address any queries or bug reports to Erik Bergstrom at ebergstr@eng.ucsd.edu
