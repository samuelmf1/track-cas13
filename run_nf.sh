#!/bin/bash
#SBATCH --job-name=nf-TRK                      # Job name
#SBATCH --mem=8G                               # Job memory request
#SBATCH --cpus-per-task=1                      # number of cpu per task
#SBATCH --time=96:00:00

module load nextflow

nextflow run main.nf -resume -profile slurm
