#!/bin/bash
#SBATCH -J group_level_pred_disease
#SBATCH -o group_level_pred_disease.o%j
#SBATCH -e group_level_pred_disease.o%j
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

# Calculate gorup level prediction accuracy
Rscript 11i_group_level_pred_disease.R

# Plot group level results
Rscript 11j_group_pred_plots_disease.R

conda deactivate
