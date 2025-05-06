#!/bin/bash
#SBATCH --job-name=protein_snp
#SBATCH --output=out/slurm_output_%j.txt
#SBATCH --error=out/slurm_error_%j.txt
#SBATCH --cpus-per-task=1
#SBATCH --mem=30G
#SBATCH --partition=jclarke,batch
#SBATCH --time=02:00:00

module load bcftools/1.17
ml anaconda
conda activate vc_env

Rscript effect_for_multi_env.R
