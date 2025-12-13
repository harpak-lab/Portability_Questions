#!/bin/bash
#SBATCH -J pc_dist_fst_plots_300K
#SBATCH -o pc_dist_fst_plots_300K.o%j
#SBATCH -e pc_dist_fst_plots_300K.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Get the proxy for the 3 1000 Genomes centroids
Rscript 15a_get_kgp_proxy_300K.R

# Make plots
Rscript 15b_pc_dist_fst_plots_300K.R

conda deactivate
