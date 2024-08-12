#!/bin/sh
#SBATCH -J ind_pred_plots
#SBATCH -o ind_pred_plots.o%j
#SBATCH -e ind_pred_plots.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

Rscript 08e_ind_pred_plots.R
