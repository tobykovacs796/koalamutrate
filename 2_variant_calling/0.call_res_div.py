# -*- coding: utf-8 -*-
"""
This script call variants in BR RESOLUTION for all individuals per chromosomes
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
ref_dir = "{}/{}/ref_fasta/{}.fasta".format(path, sp, refGenome)
bam_dir = "{}/{}/bam_files/".format(path, sp)
vcf_dir = "{}/{}/vcf_files/".format(path, sp)
#chrom_dir = "{}/{}/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}.uniq.rmdup.bam".format(path, sp, name)

# The function:
def call_var(ref, in_bam, out_dir, subset, name):
    """Happlotype caller function to call variants for each samples"""
    call_cmd = "#PBS -J {}-1906:10 \n".format(subset)
    call_cmd += "#PBS -N gVCFs_{}_{} \n".format(name, subset)
    call_cmd += "\n"
    call_cmd += "format_number() { \n"
    call_cmd += "	printf \"%04d\" \"$1\" \n"
    call_cmd += "   } \n"
    call_cmd += "index=$(format_number $PBS_ARRAY_INDEX) \n" 
    call_cmd += "\n"   
    call_cmd += "module load gatk \n"
    call_cmd += "\n"
    call_cmd += "gatk HaplotypeCaller "
    call_cmd += "-R {} ".format(ref)
    call_cmd += "-I {} ".format(in_bam)
    call_cmd += "-O {}_MSTS0100${{index}}.1_res.g.vcf ".format(out_dir)
    call_cmd += "-ERC BP_RESOLUTION "
    call_cmd += "-L MSTS0100${index}.1 "
    call_cmd += "--dont-use-soft-clipped-bases "
    call_cmd += "--native-pair-hmm-threads 1 "
    call_cmd += "--tmp-dir {}/{}/ ".format(path, sp)
    """Create a .pbs files with the calling variant functions."""
    file = open('{}_call_res_g_subset{}.pbs'.format(out_dir, subset),'w')  
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
    file.write('#PBS -l select=1:ncpus=1:mem=10GB \n')
    file.write('#PBS -l walltime=20:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write(call_cmd)
    file.close()
    ##"""Submit the .pbs to the server"""
    sub_cmd = "qsub -o {}_call_res_g_{}.out {}_call_res_g_subset{}.pbs".format(out_dir, subset, out_dir, subset)
    subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# For all merged bam files: keep uniq reads and remove duplicates.

vcf_files_dir = open("{}/{}/vcf_files.txt".format(path, sp), "w")
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}_res.g.vcf".format(vcf_dir, name)):
        print("\t The res.g.vcf file for {} already exists --> CALL VARIANT DONE".format(name))
    else:
        print("\t The res.g.vcf file for {} doesn't exist --> submit the function".format(name))
        for i in range(1, 11):
            call_var(subset=i,ref=ref_dir, in_bam=bamfile_dir[name], out_dir="{}{}".format(vcf_dir, name),name="{}".format(name))
vcf_files_dir.close()


## NOTE : range in python is not inclusive of the last value, hence if you want to include numbers 1-10 you have to say  for i in range(1, 11):
