# -*- coding: utf-8 -*-
"""
This script keeps only the uniq mapped reads,
and remove the duplicates.
It finish with the indexing of this final file.

"""
##################################################
# What you need ##################################
##################################################

# Packages:
import subprocess
import os
from variable import *

# Directories:
directory = "{}/{}/bam_files/".format(path, sp)

# Dictionary made of tuples of sample name and merged bamfile:
f = open('{}/{}/bam_files_directories.txt'.format(path, sp))
bamfile_dir = {}
for line in f:
    name = line.split()[0]
    if name not in bamfile_dir:
        bamfile_dir[name] = []
    bamfile_dir[name] = "{}/{}/bam_files/{}.bam".format(path, sp, name)


# The function:
def uniq_rmdup(input_bam, output_u, output_d, direct, sp):
    """The samtools keep uniq function."""
    uniq_cmd = "#samtools view -h {} | ".format(input_bam)
    uniq_cmd += "#grep -v -e 'XA:Z:' -e 'SA:Z:' | "
    uniq_cmd += "#samtools  view -b > {}".format(output_u)
    
    """The remove duplicate function"""
    rmdup_cmd = "samtools view -b -F 0x400 -o {} {}  ".format(output_d, output_u)
      
    """The samtools index function."""
    index_cmd = "samtools index {}".format(output_d)
    
    """Create a .pbs files with the keep uniq and remove duplicates functions."""
    file = open('{}_uniq_rmdup_idx.pbs'.format(direct),'w')
    file.write('#!/bin/bash \n')
    file.write('#PBS -P RDS-FSC-Nocticola-RW \n')
    file.write('#PBS -N remove_dup \n')
    file.write('#PBS -l select=1:ncpus=6:mem=100GB \n')
    file.write('#PBS -l walltime=20:00:00 \n')
    file.write('#PBS -q defaultQ \n')
    file.write('#PBS -M tkov3622@uni.sydney.edu.au \n')
    file.write('#PBS -m abe \n')
    file.write('module load gatk \n')
    file.write('module load samtools \n')
 #!#   file.write('module load picard \n')
    file.write(uniq_cmd)
    file.write('\n')
    file.write('echo \" Keeping unique IS DONE ##########################################################\" \n')
    file.write(rmdup_cmd)
    file.write('\n')
    file.write(index_cmd)
    file.write('\n')
    file.close()
    ##"""Submit the .pbs to the server"""


##################################################
# What you run  ##################################
##################################################

# For all merged bam files: keep uniq reads and remove duplicates.
for name in bamfile_dir:
    print(name)
    if os.path.exists("{}{}.uniq.bam".format(directory, name)):
        print("The uniq file for {} already exists".format(name))
        if os.path.exists("{}{}.uniq.rmdup.bam".format(directory, name)):
            print("The rmdup file for {} also exists ---> DONE".format(name))
        else:
            print("BUT the rmdup file for {} doesn't exists there is a problem".format(name))
    else:
        print("The uniq file for {} doesn't exist --> submit the function".format(name))
        uniq_rmdup(input_bam="{}{}.bam".format(directory, name), output_u="{}{}.uniq.bam".format(directory, name), output_d= "{}{}.uniq.rmdup.bam".format(directory, name), direct="{}{}".format(directory, name), sp=sp)
