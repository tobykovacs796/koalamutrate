# -*- coding: utf-8 -*-
"""
This script gather information about the coverage of the bam files.
"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
in_dir = "{}/{}/bam_files/".format(path, sp)
out_dir = "{}/{}/bam_files/coverage/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}.bam".format(path, sp, name)


# The function:
def coverage(in_dir, seq, out_dir):
    """The samtools function"""
    cov_cmd = "samtools depth {}{}.uniq.rmdup.bam 1> {}{}_coverage.txt".format(in_dir, seq, out_dir, seq)
    print("\t samtools depth for coverage --> {}{}_coverage.txt".format(out_dir, seq))
    """Summarize coverage"""
    sum_cmd = "python {}{}_coverage_python.py".format(out_dir, seq)
    print("\t Summarize the depth with a python script --> {}summary_coverage.txt".format(out_dir))
    """Prepare"""
    subprocess.call("cp coverage_python.py {}{}_coverage_python.py".format(out_dir, name), shell=True)
    subprocess.call("echo \"average_column(csv='{}{}_coverage.txt', name='{}')\" >> {}{}_coverage_python.py \n".format(out_dir, name, name, out_dir, name), shell=True)
    subprocess.call("sed -i '1s|^|sp=\"{}\"\\n|' {}{}_coverage_python.py \n".format(sp, out_dir, name), shell=True)
    subprocess.call("sed -i '1s|^|path=\"{}\"\\n|' {}{}_coverage_python.py \n".format(path, out_dir, name), shell=True)
    """Create a .pbs files."""
    file = open('{}{}_coverage.pbs'.format(out_dir, seq),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
    file.write('#PBS -N coverage \n')
    file.write('#PBS -l select=1:ncpus=2:mem=30GB \n')
    file.write('#PBS -l walltime=6:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write('\n')
    file.write('module load samtools \n')
    file.write('module load anaconda \n')      
    file.write('\n')
    file.write(cov_cmd)
    file.write('\n')
    file.write(sum_cmd)
    file.write('\n')
    file.close()
    """Submit the .pbs to the server"""

# The function for getting the number of reads:
def nbread_sum(seq, direct):
    """The samtools summary function."""
    nb_merged_cmd = "samtools view {}{}.bam | wc -l".format(direct, seq)
    nb_uniq_cmd = "samtools view {}{}.uniq.bam | wc -l".format(direct, seq)
    nb_rmdup_cmd = "samtools view {}{}.uniq.rmdup.bam | wc -l".format(direct, seq)
    print("\t Run number of reads in summary")
    """Create a .pbs files with the samtools summary function."""
    file = open('{}summary/{}_nb_reads.pbs'.format(direct, seq),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
    file.write('#PBS -N nbreads \n')
    file.write('#PBS -l select=1:ncpus=3:mem=50GB \n')
    file.write('#PBS -l walltime=10:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write('\n')
    file.write('module load samtools \n')       
    file.write('\n')
    file.write('SEQ={} \n'.format(seq))
    file.write('NB_MERGED=$({}) \n'.format(nb_merged_cmd))
    file.write('NB_UNIQ=$({}) \n'.format(nb_uniq_cmd))
    file.write('NB_RMDUP=$({}) \n'.format(nb_rmdup_cmd))
    file.write('echo $SEQ $NB_MERGED . $NB_UNIQ $NB_RMDUP >> {}summary/nb_reads.txt \n'.format(direct))
    file.close()

# The function for getting the name of chromosomes:
def chr_name(seq, direct):
    """The samtools summary function."""
    chr_cmd = "samtools idxstats {}{}.uniq.rmdup.bam > {}summary/chrom_name.txt".format(direct, seq, direct)
    """Create a .pbs files with the chrom name function."""
    file = open('{}summary/chrom_name.pbs'.format(direct),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
    file.write('#PBS -N chrom_names \n')
    file.write('#PBS -l select=1:ncpus=16:mem=64GB \n')
    file.write('#PBS -l walltime=10:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write('\n')
    file.write('module load samtools \n')  
    file.write('\n')     
    file.write(chr_cmd)
    file.write('\n')
    file.close()

##################################################
# What you run  ##################################
##################################################
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}.uniq.rmdup.bam".format(in_dir, name)):
        print("\t The uniq_rmdup file for {} exists".format(name))
        # Count the number of reads:
        nbread_sum(seq=name, direct=in_dir)
    if os.path.exists("{}{}_coverage.txt".format(out_dir, name)):
        print("\t The coverage of {} has already been calculated --> NOT CALCULATED".format(name))
    else:
        print("\t The coverage of {} has NOT been calculated".format(name))
        coverage(in_dir=in_dir, seq=name, out_dir=out_dir)

chr_name(seq=name, direct=in_dir)
print("\n Look at the name of chromosome only for {} sequence".format(name))
