# -*- coding: utf-8 -*-
"""
This script write calculate the fasle negative rate:
       - prepare the file
       - write the R script to find the number of true variant outside the threshold

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
direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
direct_handling="{}/{}/vcf_handling/".format(path, sp)
direct = "{}/{}/".format(path, sp)

# Dictionary:
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

# Import depth
depth_mean = pd.read_csv('{}depth.txt'.format(direct),sep=' ', index_col=None, header=None)

# The functions:
def fnr(direct, direct_handling, direct_denovo, name, GQ_lim, DP_min, DP_max, AB_max, AB_min, father, mother, offspring):
	"""Write the start of the Rscript"""
	file = open('{}false_neg_rate_{}.r'.format(direct_denovo, name),'w')
	file.write('# Sample: \n')
	file.write('name=\"{}\" \n'.format(name))
	file.write('\n')
	file.write('# Directories: \n')
	file.write('direct_denovo=\"{}\" \n'.format(direct_denovo))
	file.write('direct=\"{}\" \n'.format(direct))
	file.write('\n')
	file.write('\n')
	file.write('# Filters: \n')
	file.write('AB_max={} \n'.format(AB_max))
	file.write('AB_min={} \n'.format(AB_min))
	file.write('\n')
	file.close()
	"""Catenate with the rest"""
	cat_cmd = "cat {}rscript_false_neg_rate.R >>{}false_neg_rate_{}.r".format(direct_denovo, direct_denovo, name)
	"""BCFtools"""
	
	bcf1_cmd = "bcftools view -i '(FMT/GT[{}]=\"AA\" && FMT/GT[{}]=\"RR\" && FMT/GT[{}]=\"het\")||(FMT/GT[{}]=\"RR\" && FMT/GT[{}]=\"AA\" && FMT/GT[{}]=\"het\")' {}inter_vcf/genotype_genomicDBI_{}_snp_filt.g.vcf >> {}neg_rate_1_{}.g.vcf".format(father,mother,offspring,father,mother,offspring,direct_handling, name, direct_denovo, name)
	bcf2_cmd= "bcftools view -i 'MIN(FMT/DP)>={} & MIN(FMT/GQ)>={} & MAX(FMT/DP)<={}' {}neg_rate_1_{}.g.vcf >> {}neg_rate_2_{}.g.vcf".format(DP_min, GQ_lim, DP_max, direct_denovo, name, direct_denovo, name)
	"""Convert to table"""
	vtt_cmd = "gatk VariantsToTable "
	vtt_cmd += "-V {}neg_rate_2_{}.g.vcf ".format(direct_denovo, name)
	vtt_cmd += "-F CHROM -F POS -F TYPE -F REF -F ALT -F FILTER -GF GT -GF AD -GF DP -GF GQ -GF PL -GF SAC "
	vtt_cmd += "-O {}output.table.fnr.{} ".format(direct_denovo, name)
	vtt_cmd += "--show-filtered"
	"""Run r script"""
	r_cmd = "Rscript {}false_neg_rate_{}.r \n".format(direct_denovo, name)
	"""Create a .pbs files with the false negative rate functions."""
	file = open('{}{}_fnr.pbs'.format(direct_denovo, name),'w')
	file.write('#!/bin/bash \n')
	file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
	file.write('#PBS -l select=1:ncpus=1:mem=64GB \n')
	file.write('#PBS -l walltime=1:00:00 \n')
	file.write('#PBS -q defaultQ \n')
	file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
	file.write('#PBS -m abe \n')
	file.write('\n')
	file.write('module load bcftools/1.9 \n')
	file.write('module load R \n')
	file.write('module load gatk \n')
	file.write(bcf1_cmd)
	file.write('\n')
	file.write(bcf2_cmd)
	file.write('\n')
	file.write('## Convert to table: \n')
	file.write(vtt_cmd)
	file.write('\n')
	file.write('## Catenate the Rscript for samples and remaining part: \n')
	file.write(cat_cmd)
	file.write('\n')
	file.write('## Run the Rscript: \n')
	file.write(r_cmd)
	file.write('\n')
	file.close()
	##"""Submit the .pbs to the server"""
	sub_cmd = "qsub {}{}_fnr.pbs".format(direct_denovo, name, direct_denovo, name)
	subprocess.call(sub_cmd, shell=True)


##################################################
# What you run  ##################################
##################################################

# Find the false negative rate per trio:
for name, trios in trio_dir.items():
    for trio in trios:
        offspring, father, mother = trio  # Unpack the trio
        
        # Now you can use offspring, father, and mother in your calculations
        mean_dp = depth_mean.loc[depth_mean[0] == offspring][1]
        DP_min_calc = int(round(eval(DP_min)))
        DP_max_calc = int(round(eval(DP_max)))
        fnr(direct=direct, direct_denovo=direct_denovo, direct_handling=direct_handling, name=name, GQ_lim=GQ_lim, DP_min=DP_min_calc, DP_max=DP_max_calc, AB_max=AB_max, AB_min=AB_min, father=father, mother=mother, offspring=offspring)
        print("Finding the false negative rate for trio {}".format(name))

