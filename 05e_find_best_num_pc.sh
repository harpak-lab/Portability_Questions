#!/bin/bash
#SBATCH -J find_best_num_pc
#SBATCH -o find_best_num_pc.o%j
#SBATCH -e find_best_num_pc.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Format Fst files
Rscript 05f_final_fst_formatting.R

# Identify the best number of PCs to use
Rscript 05g_best_num_pc.R

conda deactivate
