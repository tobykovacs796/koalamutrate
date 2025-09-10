# -*- coding: utf-8 -*-
"""
This script catenate the GVCF files per chromosome for all individuals.
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


# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}.uniq.rmdup.bam".format(path, sp, name)

# And VCF
f = open('{}/{}/vcf_files.txt'.format(path, sp))
vcf_dir = {}
for line in f:
    name = line.split()[0]
    if name not in vcf_dir:
        vcf_dir[name] = []
    vcf_dir[name] = "{}/{}/vcf_files/{}".format(path, sp, name)


# Import chrom names:
chrom_name = pd.read_csv('{}chromosomes.txt'.format(chrom_dir),sep=' ', index_col=None, header=None)

# The function:
def combine(direct):
    """Combine all samples with GenomicsDBImport"""        
    combine_cmd = "#PBS -N comb_gVCFs_subset10 \n"
    combine_cmd += "#PBS -J 1-1906:10 \n"
    combine_cmd += "\n"
    combine_cmd += "module load gatk"
    combine_cmd += "\n"
    combine_cmd += "format_number() { \n"
    combine_cmd += "	printf \"%04d\" \"$1\" \n"
    combine_cmd += "   } \n"
    combine_cmd += "\n"
    combine_cmd += "index=$(format_number $PBS_ARRAY_INDEX) \n" 
    combine_cmd += "\n"
    combine_cmd += "gatk GenomicsDBImport " 
    combine_cmd += "--variant {}TarZoo_M_B70170_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_M_46879_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}TarZoo_F_B80238_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_F_46850_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_M_46851_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_F_46861_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}TarZoo_F_B20093_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}TarZoo_M_B60339_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_M_46878_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_F_46860_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_F_46865_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--variant {}Featherdale_F_46857_MSTS0100${{index}}.1_res.g.vcf ".format(direct)
    combine_cmd += "--tmp-dir {}/{}/ ".format(path, sp)
    combine_cmd += "--genomicsdb-workspace-path {}/genomicDBI_MSTS0100${{index}}.1 ".format(direct)
    combine_cmd += "-L MSTS0100${index}.1 "
    """Create a .pbs files with the combine variant functions."""
    file = open('{}combine_genomicDBImport_subset10.pbs'.format(direct),'w')
    
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-BacCock-RW \n')
    file.write('#PBS -l select=1:ncpus=1:mem=20GB \n')
    file.write('#PBS -l walltime=2:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write(combine_cmd)
    file.write('\n')
    file.close()


