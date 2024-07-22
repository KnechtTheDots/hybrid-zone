#!/bin/bash
#SBATCH --job-name=<j<ob name here>>
#SBATCH --output=snakemake_out/<output file here>
#SBATCH --time=2-00:00:00
#SBATCH --mail-user=jknecht1@binghamton.edu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1

eval "$(conda shell.bash hook)"

conda activate hyb-zone

snakemake --profile slurm
