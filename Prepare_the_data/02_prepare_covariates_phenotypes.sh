#!/bin/bash
#SBATCH -J prepare_covariates_phenotypes
#SBATCH -o prepare_covariates_phenotypes.o%j
#SBATCH -e prepare_covariates_phenotypes.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Creates covariate file with age, sex, age*sex, age^2, age^2 * sex and PC1, ..., PC20
Rscript 02a_prepare_covariates.R

# Creates a file contining all 15 phenotypes
Rscript 02b_prepare_phenotypes.R

conda deactivate
