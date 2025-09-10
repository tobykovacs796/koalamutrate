# -*- coding: utf-8 -*-
"""
This script calls variants with mpilup for each 
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
import pandas as pd
from variable import *


# Import species names:
#species = pd.read_csv('sp_here.txt', sep=' ', index_col=None, header=None)
#pedigree = pd.read_csv('pedigree.txt', sep='\t', index_col=None, header=None)
#genome = pd.read_csv('ref.txt', sep='\t', index_col=None, header=None)

##################################################
# What you run  ##################################
##################################################

# Directories:
direct_denovo="{}/{}/de_novo_mutation/".format(path, sp)
direct_handling="{}/{}/vcf_handling/".format(path, sp)
direct_bam="{}/{}/bam_files/".format(path, sp)
direct = "{}/{}/".format(path, sp)
ref_dir = "{}/{}/ref_fasta/{}.fasta".format(path, sp, refGenome)

    # Dictionary:
pedigree_sp=pd.read_csv('{}/pedigree.ped'.format(direct), sep=' ', index_col=None, header=None)

for sample in range(1,len(pedigree_sp)):
	fa = pedigree_sp.iloc[sample,2]
	mo = pedigree_sp.iloc[sample,3]
	off = pedigree_sp.iloc[sample,1]
	denovo_to_check=pd.read_csv('{}data_denovo_{}.tab'.format(direct_denovo, off), sep='\t', index_col=None)
	file = open('{}/{}_samtools.pbs'.format(direct_denovo, off),'w')
	file.write('#!/bin/bash \n')
	file.write('#PBS -P RDS-FSC-BacCock-RW \n')
	file.write('#PBS -l select=1:ncpus=1:mem=19GB \n')
	file.write('#PBS -l walltime=1:00:00 \n')
	file.write('#PBS -q defaultQ \n')
	file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
	file.write('#PBS -m abe \n')
	file.write('\n')  
	file.write('module load samtools/1.2 \n')
	file.write('module load bcftools/1.2 \n')
	for line in range(0,len(denovo_to_check)):
		chrom=denovo_to_check.iloc[line,0]
		pos=denovo_to_check.iloc[line,1]
		file.write('samtools mpileup -ugf {} -r {}:{}-{} {}{}.uniq.rmdup.bam {}{}.uniq.rmdup.bam {}{}.uniq.rmdup.bam | bcftools call -m | tail -1 | cut -f1,2,4,5,10,11,12 >> {}{}_samtools.txt \n'.format(ref_dir, chrom, pos, pos, direct_bam, fa, direct_bam, mo, direct_bam, off, direct_denovo, off))
	file.close()


