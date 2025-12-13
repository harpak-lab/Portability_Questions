#!/bin/sh
#SBATCH -J calc_ma_counts
#SBATCH -o calc_ma_counts.o%j
#SBATCH -e calc_ma_counts.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

Rscript 08c_calc_ma_counts.R
