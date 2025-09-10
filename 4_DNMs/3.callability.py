# -*- coding: utf-8 -*-
"""
This script generates the callability
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
direct="{}/{}/".format(path, sp)

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

##################################################
# What you run  ##################################
##################################################

    # Create a file with the callability and import:
for name in trio_dir:
	call_cmd= "echo \"{} $(awk 'NR==4' {}callability_{}.out)\" >> {}callability.txt".format(name, direct_denovo, name, direct_denovo)
	subprocess.call(call_cmd, shell=True)
