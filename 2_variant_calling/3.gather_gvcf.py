# -*- coding: utf-8 -*-
"""
This script gather:
    the genotypes per chromosomes into one files
    the back combine per chromosomes into one file
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *
import pandas as pd

# Directories:
direct = "{}/{}/vcf_files/".format(path, sp)
chrom_dir = "{}/{}/".format(path, sp)

# Import chrom names:
chrom_name = pd.read_csv('{}chromosomes.txt'.format(chrom_dir),sep=' ', index_col=None, header=None)

# List of the files:
back_dir = []
geno_dir = []
for line in range(0,nb_chrom):
    chrom=chrom_name.loc[line,1]
    back_dir.append("back_combine_genomicDBI_{}.g.vcf".format(chrom))
    geno_dir.append("genotype_genomicDBI_{}.g.vcf".format(chrom))


# The function:
def gather(direct, list_file, output, what):
    """Gather each chrom together"""
    gather_cmd = "#PBS -N gather_gVCFs \n"
    gather_cmd += "\n"
    gather_cmd += "module load gatk"
    gather_cmd += "\n"
    gather_cmd += "gatk GatherVcfs "
    for i in list_file:
        gather_cmd += "-I {}{} ".format(direct, i)
    gather_cmd += "-O {} ".format(output)
    """Create a .pbs files with the gather functions."""
    file = open('{}gather_genomicDBImport_{}.pbs'.format(direct, what),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-BacCock-RW \n')
    file.write('#PBS -l select=1:ncpus=1:mem=124GB \n')
    file.write('#PBS -l walltime=12:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write(gather_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .pbs to the server"""


##################################################
# What you run  ##################################
##################################################

# Gather if everything exist for all chromosomes:
list_exist=[]
for line in range(0, nb_chrom):
    chrom=chrom_name.loc[line,1]
    list_exist.append(os.path.exists("{}back_combine_genomicDBI_{}.g.vcf".format(direct, chrom)))
if all(list_exist):
    print("\t All back_combine directories exist --> back combine done")
    mv_bc= "mv {}back_combine_genomicDBImport_* {}back_com.log/".format(direct, direct)
    subprocess.call(mv_bc, shell=True)
    print("\t Move the back combine log files")
    gather(direct=direct, list_file=back_dir, output="{}back_combine_genomicDBI_gather.g.vcf".format(direct), what="back_combine")
    print("Gather back combine")

list_exist=[]
for line in range(0, nb_chrom):
    chrom=chrom_name.loc[line,1]
    list_exist.append(os.path.exists("{}genotype_genomicDBI_{}.g.vcf".format(direct, chrom)))
if all(list_exist):
    print("\t All genotypes directories exist --> genotype done")
    mv_geno= "mv {}genotype_genomicDBImport_* {}back_com.log/".format(direct, direct)
    subprocess.call(mv_geno, shell=True)
    print("\t Move the genotype log files")
    gather(direct=direct, list_file=geno_dir, output="{}genotype_genomicDBI_gather.g.vcf".format(direct), what="genotype")
    print("Gather genotypes")
