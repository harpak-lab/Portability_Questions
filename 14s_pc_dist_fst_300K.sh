#!/bin/bash
#SBATCH -J pc_dist_fst_300K
#SBATCH -o pc_dist_fst_300K.o%j
#SBATCH -e pc_dist_fst_300K.o%j
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

# Calculate weighted PC distance
Rscript 14t_pc_dist_calc_300K.R

# Edit the fam file for Fst calculations
Rscript 14u_edit_merged_fam_300K.R

# Write scripts for Fst calculations
python 14v_write_fst_scripts_300K.py

conda deactivate

# Then run the following in the command lines
# cd temp_fst_300K_path
# for i in $(seq 0 49);
# do
#   sbatch run_$i.sh
# done
