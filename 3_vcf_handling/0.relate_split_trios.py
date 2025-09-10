# -*- coding: utf-8 -*-
"""
This script test the relatedness between individuals and create a genotype/back_combine files per trios.
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
vcf_dir="{}/{}/vcf_files/".format(path, sp)
direct = "{}/{}/vcf_handling/".format(path, sp)
ref="{}/{}/ref_fasta/{}.fasta".format(path, sp,refGenome)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/pedigree.ped'.format(path, sp))
trio_dir = {}
for line in f:
    off = line.split()[1]
    fa = line.split()[2]
    mo = line.split()[3]
    name = off
    if name not in trio_dir:
        trio_dir[name] = []
    trio_dir[name].append((off, fa, mo))



# The function:
def relatedness(vcf_in, output, direct):
	"""Relatedness between individuals"""
	rn_cmd = "vcftools --vcf {} --out {} --relatedness2".format(vcf_in, output)
	"""Create a .pbs files with the relatedness functions."""
	file = open('{}relatedness.pbs'.format(direct),'w')
	file.write('#!/bin/bash \n')
	file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
	file.write('#PBS -l select=1:ncpus=1:mem=56GB \n')
	file.write('#PBS -l walltime=1:00:00 \n')
	file.write('#PBS -q defaultQ \n')
	file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
	file.write('#PBS -m abe \n')
	file.write('\n')
	file.write('module load vcftools')
	file.write('\n')
	file.write('\n')
	file.write(rn_cmd)
	file.write('\n')
	file.close()
	##"""Submit the .pbs to the server"""


def select_trio(ref, geno_file, off, fa, mo, direct, output_trio, cpu, mem_j, mem_r, time):
	"""Select a trio"""
	sel_cmd = "gatk SelectVariants "
	sel_cmd += "-R {} ".format(ref)
	sel_cmd += "-V {} ".format(geno_file)
	sel_cmd += "-sn {} ".format(off)
	sel_cmd += "-sn {} ".format(fa)
	sel_cmd += "-sn {} ".format(mo)
	sel_cmd += "-O {}{}.g.vcf ".format(direct, output_trio)
	"""Create a .pbs files with the select trio functions."""
	file = open('{}{}.pbs'.format(direct, output_trio),'w')
	file.write('#!/bin/bash \n')
	file.write('#PBS -P RDS-FSC-BacCock-RW \n')
	file.write('#PBS -l select=1:ncpus=1:mem=20GB \n')
	file.write('#PBS -l walltime=1:00:00 \n')
	file.write('#PBS -q defaultQ \n')
	file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
	file.write('#PBS -m abe \n')
	file.write('\n')
	file.write('module load gatk')
	file.write('\n')
	file.write(sel_cmd)
	file.write('\n')
	file.close()
    ##"""Submit the .pbs to the server"""


##################################################
# What you run  ##################################
##################################################

##Test for relatedness:
relatedness(vcf_in="{}genotype_genomicDBI_gather.g.vcf".format(vcf_dir), output="{}{}".format(direct, sp), direct=direct)
print("Test for relatedness between all individuals \n")

# Create per trio files:
for name in trio_dir:
    off = trio_dir[name][0][0]
    fa = trio_dir[name][0][1]
    mo = trio_dir[name][0][2]
    select_trio(ref=ref, geno_file="{}genotype_genomicDBI_gather.g.vcf".format(vcf_dir), off= off, fa=fa, mo=mo, direct=direct, output_trio="genotype_genomicDBI_{}".format(name), cpu=1, mem_j=60, mem_r=64, time=10)
    print("VCF file created for trio {} --> genotype_genomicDBI_{}.g.vcf".format(name, name))
    select_trio(ref=ref, geno_file="{}back_combine_genomicDBI_gather.g.vcf".format(vcf_dir), off= off, fa=fa, mo=mo, direct=direct, output_trio="back_combine_genomicDBI_{}".format(name), cpu=1, mem_j=90, mem_r=100, time=20)
    print("VCF file created for trio {} --> back_combine_genomicDBI_{}.g.vcf".format(name, name))

print("{} files created for genotype_genomicDBI".format(len(trio_dir)))
print("{} files are created for back_combine_genomicDBI".format(len(trio_dir)))
