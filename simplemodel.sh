#!/bin/bash

# setting the name of the job
#SBATCH --job-name=simplemodel

# setting home directory
#SBATCH -D /home/astock/chileanatn/sims

# setting standard error output
#SBATCH -e /home/astock/chileanatn/sims/slurm_log/sterror_%j.txt

# setting standard output
#SBATCH -o /home/astock/chileanatn/sims/slurm_log/stdoutput_%j.txt

# setting account
#SBATCH -A fvaldovigrp

# setting medium priority
#SBATCH -p med2

# setting the maximum time
#SBATCH -t 12:00:00

# set maximum memory
#SBATCH --mem=100GB

# mail alerts at beginning and end of job
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END


# send mail here
#SBATCH --mail-user=astockdill@ucdavis.edu

# now we'll print out the contents of the R script to the standard output file
# cat simulations.jl
# echo "ok now for the actual standard output"

# load julia
module load julia

echo "finished loading, starting file"
# run script
#julia -p $SLURM_NTASKS simulations.jl

# Read parameters
#IFS=' ' read -r guild B0 TI <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" params_list.txt)"
#echo "Running: Guild=$guild, B0=$B0, TI=$TI"

#echo $guild $B0 $TI $T

julia --project=. simplemodel.jl

