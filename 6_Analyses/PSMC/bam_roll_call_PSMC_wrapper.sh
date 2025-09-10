#!/bin/bash

# Check if the argument is provided
if [ $# -ne 1 ]; then
    echo "Please provide an argument."
    echo "Usage: $0 <argument>"
    exit 1
fi

# Create or overwrite the output file
> "bams_$1.txt"

# Loop through each BAM file in the directory
for bam_file in *.bam; do
    # Extract the individual's name from the BAM file name
    individual=$(basename "$bam_file" .bam)
    # Append the individual's name to the output file
    echo "$individual" >> "bams_$1.txt"
done

# Count the number of individuals in the bams file
num_individuals=$(wc -l < "bams_$1.txt")

# Generate the script containing the specified content
script_file="2.psmc_$1.pbs"

cat > "$script_file" <<EOF
#!/bin/bash
#PBS -P RDS-FSC-Nocticola-RW
#PBS -N PSMC_filtering
#PBS -l select=1:ncpus=4:mem=80GB
#PBS -l walltime=100:00:00
#PBS -q defaultQ
#PBS -M tkov3622@uni.sydney.edu.au
#PBS -m abe
#PBS -J 1-$num_individuals

params=\$(sed "\${PBS_ARRAY_INDEX}q;d" /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/bams/$1/bams_$1.txt)
param_array=( \$params )

# Load modules
module load samtools
module load bcftools
module load psmc
module load gnuplot

# Set working directory
cd /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/bams/$1

# Run it
samtools view -@ 4 -bh -F 0x400 -L /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/ref_genomes/koala.50k_noX.bed -o \${param_array[0]}.50k_noX.bam \${param_array[0]}.bam
samtools index -@ 4 \${param_array[0]}.50k_noX.bam
samtools flagstat -@ 4 \${param_array[0]}.50k_noX.bam > \${param_array[0]}.flagstat
samtools depth \${param_array[0]}.50k_noX.bam | awk '{sum+=\$3} END { print "Average = ",sum/NR}' > \${param_array[0]}.50k_noX.coverage.txt

cov=\$(awk '{printf("%d"), \$3}' \${param_array[0]}.50k_noX.coverage.txt)
mincov=\$((\$cov / 3))
maxcov=\$((\$cov * 2))
samtools mpileup -A -C50 -uf /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/ref_genomes/GCA_002099425.1_phaCin_unsw_v4.1_genomic.fasta \${param_array[0]}.50k_noX.bam | bcftools call -c - |
vcfutils.pl vcf2fq -d \${mincov} -D \${maxcov} | gzip > \${param_array[0]}.con.fq.gz

fq2psmcfa -q20 \${param_array[0]}.con.fq.gz > \${param_array[0]}.psmcfa

mkdir /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/psmc/$1
cp \${param_array[0]}.psmcfa /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/psmc/$1
cd /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/PSMC/psmc/$1

## run the PSMC which is almost instantaneous

psmc -N25 -t15 -r5 -p "4+25*2+4+6" -o \${param_array[0]}.4.psmc \${param_array[0]}.psmcfa
psmc_plot.pl -R -g 7 -u 0.00000000612 \${param_array[0]}.612.4 \${param_array[0]}.4.psmc


psmc -N25 -t15 -r5 -p "2+2+25*2+4+6" -o \${param_array[0]}.22.psmc \${param_array[0]}.psmcfa
psmc_plot.pl -R -g 7 -u 0.00000000612 \${param_array[0]}.612.22 \${param_array[0]}.22.psmc

psmc -N25 -t15 -r5 -p "1+1+1+1+25*2+4+6" -o \${param_array[0]}.1111.psmc \${param_array[0]}.psmcfa
psmc_plot.pl -R -g 7 -u 0.00000000612 \${param_array[0]}.612.1111 \${param_array[0]}.1111.psmc


EOF

# Make the generated script file executable
chmod +x "$script_file"
