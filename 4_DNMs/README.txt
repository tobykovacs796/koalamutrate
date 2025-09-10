Set up the right variable in variable.py

Calculate callability, number of de novo mutation, false negative rate and mutation rate:
    python 0.callability.py --> give the number of callable sites (require bcftools version 1.9)
    python 1.negative_rate.py --> give the number of true heterozygotes outside the AB boundaries (require bcftools version 1.9)
    python 2.nb_de_novo.py --> give the number of mutation per trio
    python 3.callability.py --> creates the callability file required in next step
    python 4.rate_samtools.py --> for all samples, produce variant calling with samtools mpileup and bcftools call (samtools and bcftools version 1.2)
    python 5.rate_samtools_denovo.py --> re calculate the number of DNMs (without the appearant FPs)
    python 6.rate_samtools_rate.py --> calculate a mutation rate as: DNM_corrected/(2xCGx(1-FNR))

