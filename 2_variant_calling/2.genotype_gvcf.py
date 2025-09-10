# -*- coding: utf-8 -*-
"""
This script genotypes and constructs back combined files
from the GenomicDBI folders VCF files per chromosome for all individuals.
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
ref = "{}/{}/ref_fasta/{}.fasta".format(path, sp, refGenome)
chrom_dir = "{}/{}/".format(path, sp)

# Import chrom names:
chrom_name = pd.read_csv('{}chromosomes.txt'.format(chrom_dir),sep=' ', index_col=None, header=None)

# The function:
def genotype(ref, sp, subset, direct):
    """Genotypes all samples from GenomicsDBImport"""
    gt_cmd = "#PBS -N comb_gVCFs_subset{} \n".format(subset)
    gt_cmd += "#PBS -J {}-1906:10 \n".format(subset)
    gt_cmd += "\n"
    gt_cmd += "module load gatk".format(subset)
    gt_cmd += "\n"
    gt_cmd += "format_number() { \n"
    gt_cmd += "	printf \"%04d\" \"$1\" \n"
    gt_cmd += "   } \n"
    gt_cmd += "\n"
    gt_cmd += "index=$(format_number $PBS_ARRAY_INDEX) \n" 
    gt_cmd += "\n"
    gt_cmd += "gatk GenotypeGVCFs "
    gt_cmd += "--tmp-dir /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/Phascolarctos_cinereus "
    gt_cmd += "-R {} ".format(ref)
    gt_cmd += "-V gendb://{}genomicDBI_MSTS0100${{index}}.1 ".format(direct)
    gt_cmd += "-G StandardAnnotation -new-qual "
    gt_cmd += "-O {}genotype_genomicDBI_MSTS0100${{index}}.1.g.vcf ".format(direct)
    """Create a .pbs files with the genotype functions."""
    file = open('{}genotype_genomicDBImport_subset{}.1.pbs'.format(direct,subset),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-BacCock-RW \n')
    file.write('#PBS -l select=1:ncpus=1:mem=70GB \n')
    file.write('#PBS -l walltime=5:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write(gt_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .pbs to the server"""
    sub_cmd = "qsub -o {}genotype_genomicDBImport_subset{}.1.out {}genotype_genomicDBImport_subset{}.1.pbs".format(direct, subset, direct, subset)
    subprocess.call(sub_cmd, shell=True)


def back_combine(ref, sp, subset, direct):
    """Combine readable from GenomicsDBImport"""
    
    bc_cmd = "#PBS -N comb_gVCFs_subset{} \n".format(subset)
    bc_cmd += "#PBS -J {}-1906:10 \n".format(subset)
    bc_cmd += "\n"
    bc_cmd += "module load gatk".format(subset)
    bc_cmd += "\n"
    bc_cmd += "format_number() { \n"
    bc_cmd += "	printf \"%04d\" \"$1\" \n"
    bc_cmd += "   } \n"
    bc_cmd += "\n"
    bc_cmd += "index=$(format_number $PBS_ARRAY_INDEX) \n" 
    bc_cmd += "\n"
    bc_cmd += "gatk SelectVariants "
    bc_cmd += "--tmp-dir /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/Phascolarctos_cinereus "
    bc_cmd += "-R {} ".format(ref)
    bc_cmd += "-V gendb://{}/genomicDBI_MSTS0100${{index}}.1 ".format(direct)
    bc_cmd += "-O {}back_combine_genomicDBI_MSTS0100${{index}}.1.g.vcf ".format(direct)
    """Create a .pbs files with the genotype functions."""
    file = open('{}back_combine_genomicDBImport_subset{}.pbs'.format(direct, subset),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-BacCock-RW \n')
    file.write('#PBS -l select=1:ncpus=1:mem=160GB \n')
    file.write('#PBS -l walltime=10:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write(bc_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .pbs to the server"""
    sub_cmd = "qsub -o {}back_combine_genomicDBImport_subset{}.1.out {}back_combine_genomicDBImport_subset{}.pbs".format(direct, subset, direct, subset)
    subprocess.call(sub_cmd, shell=True)

##################################################
# What you run  ##################################
##################################################

# For each chromosome one function:
list_exist=[]
for line in range(0, nb_chrom):
    chrom=chrom_name.loc[line,1]
    list_exist.append(os.path.exists("{}genomicDBI_{}".format(direct, chrom)))
if all(list_exist):
    print("\t All the genomicDBI directoris exist --> combine variant done")
    mv_com = "mv {}combine_genomicDBImport_* {}combine.log".format(direct, direct)
    subprocess.call(mv_com, shell=True)
    print("\t Move the combine log files")
    mv_vcf = "mv {}*_res.g.vcf* {}inter_vcf".format(direct, direct)
    subprocess.call(mv_vcf, shell=True)
    print("\t Move the inter vcf files")
    for i in range(1, 11):
        genotype(ref=ref, sp=sp, subset=i, direct=direct)
        print("Genotyping for subset{}".format(i))
        back_combine(ref=ref, sp=sp, subset=i, direct=direct)
        print("Back combining to have combine for subset{} readable".format(i))

else:
	print("Some genomicDBI directoris are missing --> PROBLEM")