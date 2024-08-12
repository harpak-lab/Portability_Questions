#!/bin/bash
#SBATCH -J compare_effect_sizes_heterozygosity_var_pgs_plots
#SBATCH -o compare_effect_sizes_heterozygosity_var_pgs_plots.o%j
#SBATCH -e compare_effect_sizes_heterozygosity_var_pgs_plots.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 2:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Make the plots
Rscript 09k_compare_effect_sizes_heterozygosity_var_pgs_plots.R

conda deactivate
