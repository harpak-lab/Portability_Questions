#!/bin/bash
#SBATCH -J townsend_income_sex_age_country
#SBATCH -o townsend_income_sex_age_country.o%j
#SBATCH -e townsend_income_sex_age_country.o%j
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

# Make plots
Rscript 12a_townsend_income_sex_age_country_plots.R

conda deactivate
