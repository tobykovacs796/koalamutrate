Set up the right variable in variable.py

Calling variant for each individual per chromosomes, combine all individuals, joint genotype and merge the chromosomes files:
    python 0.call_res_div.py --> call variable in bp resolution mode for each samples per chromosomes
    python 1.combine_gvcf_div.py --> combine all the samples together per chromosomes
    python 2.genotype_gvcf.py --> genotype per chromosomes and back combine per chromosome
    python 3.gather_gvcf.py --> to have a unique genotype and back combine for all chromosomes
    python 4.final_check.py --> check all the files are present and move intermediate files
