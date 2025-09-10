#!/bin/bash

# Define populations
populations=("N_QLD" "SEQLD_NNSW" "M_NSW" "S_NSW" "VIC_noSA")

# Loop over all unique pairs
for ((i=0; i<${#populations[@]}; i++)); do
    for ((j=i+1; j<${#populations[@]}; j++)); do
        pop1=${populations[i]}
        pop2=${populations[j]}
        log_name="log2.est_splt_${pop1}_${pop2}_50knots"
        job_name="2.est_split_${pop1}_${pop2}_50knots"
        output_dir="${pop1}_${pop2}_50knots"
        pbs_file="${job_name}.pbs"

        cat > "$pbs_file" <<EOF
#!/bin/bash

#PBS -P RDS-FSC-Nocticola-RW
#PBS -N $log_name
#PBS -l select=1:ncpus=16:mem=50GB
#PBS -l walltime=10:00:00
#PBS -m abe
#PBS -M tkov3622@hpc.sydney.edu.au
#PBS -q defaultQ

module load python/3.9.15 mpfr/4.0.2 gsl/1.16 gcc/9.1.0

cd /scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/SMCpp/full/smc/split

# Set up temporary directory
export TMPDIR=/scratch/RDS-FSC-awggdata-RW/koala_dnazoo/wgs/SMCpp/smc/tmp/job_\${PBS_JOBID}
mkdir -p "\$TMPDIR"

source smc-env/bin/activate

# Run split using estimates for each population
smc++ split --timepoints 1e3 1e6 -o $output_dir/ \\
    ../estimate/$pop1/tpe3e6_50knots/model.final.json \\
    ../estimate/$pop2/tpe3e6_50knots/model.final.json \\
    ../estimate/$pop1/data/${pop1}_*.smc.gz \\
    ../estimate/$pop2/data/${pop2}_*.smc.gz \\
    data/${pop1}_${pop2}_*.smc.gz \\
    data/${pop2}_${pop1}_*.smc.gz

# Plot all of the models
smc++ plot -c split_${pop1}_${pop2}_50knots.pdf \\
    $output_dir/model.final.json --csv -g 7
EOF

        echo "Generated PBS script: $pbs_file"
    done
done
