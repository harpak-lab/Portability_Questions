#!/bin/bash
#SBATCH -J calc_heterozygosity
#SBATCH -o calc_heterozygosity.o%j
#SBATCH -e calc_heterozygosity.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

# Prepare IDs in each bin
Rscript 08i_prepare_heterozygosity.R

# Calculate heterozygosity in each bin
for i in $(seq 1 500);
do
  $plink \
    --bfile ${scratch}/ukb_merged/ukb_merged \
    --keep data/ukb_populations/bins/bin_${i}_id.txt \
    --extract data/pgs/all_pgs_snps.txt \
    --freqx \
    --out data/heterozygosity/heterozygosity_bin_${i}
done

# Stratify the SNPs by effect sizes and calculate the mean heterozygosity
Rscript 08j_post_heterozygosity_calc.R

conda deactivate
