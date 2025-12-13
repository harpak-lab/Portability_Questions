#!/bin/bash
#SBATCH -J pc_dist_fst
#SBATCH -o pc_dist_fst.o%j
#SBATCH -e pc_dist_fst.o%j
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
Rscript 05b_pc_dist_calc.R

# Edit the fam file for Fst calculations
Rscript 05c_edit_merged_fam.R

# Write scripts for Fst calculations
python 05d_write_fst_scripts.py

conda deactivate

# Then run the following in the command lines
# cd temp_fst_path
# for i in $(seq 0 49);
# do
#   sbatch run_$i.sh
# done
