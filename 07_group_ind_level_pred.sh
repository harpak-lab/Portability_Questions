#!/bin/bash
#SBATCH -J group_ind_level_pred
#SBATCH -o group_ind_level_pred.o%j
#SBATCH -e group_ind_level_pred.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Calculate gorup level prediction accuracy
Rscript 07a_group_level_pred.R

# Calculate individual level prediction error
Rscript 07b_ind_level_pred.R

# Plot group and individual level results
Rscript 07c_group_ind_pred_plots.R

conda deactivate
