#!/bin/bash
#SBATCH -J PCA_close_far
#SBATCH -o PCA_close_far.o%j
#SBATCH -e PCA_close_far.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 2:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e 

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

$plink2 \
  --memory 70000 \
  --bfile ${scratch}/ukb_merged/ukb_merged_snps_in_pca \
  --keep data/ukb_populations/close_id.txt \
  --pca approx 20 \
  --out data/pca/close_pca

$plink2 \
  --memory 70000 \
  --bfile ${scratch}/ukb_merged/ukb_merged_snps_in_pca \
  --keep data/ukb_populations/far_id.txt \
  --pca approx 20 \
  --out data/pca/far_pca
