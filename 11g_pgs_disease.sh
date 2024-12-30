#!/bin/bash
#SBATCH -J PGS_disease
#SBATCH -o PGS_disease.o%j
#SBATCH -e PGS_disease.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e 

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

for phenotype in Alzheimer T2D Asthma
do
  for threshold in 0 1 2 3 4
  do
    # Score each individual using the SNPs below a p-value threshold and GWAS betas
    $plink2 \
      --bfile ${scratch}/pgs/CT_disease_merged \
      --extract data/pgs/${phenotype}_disease_threshold_${threshold}.txt \
      --score data/gwas_results/${phenotype}_disease_combined.glm.logistic.hybrid 3 6 14 header no-mean-imputation \
      --memory 70000 \
      --out data/pgs/${phenotype}_disease_${threshold}_scores

    $plink2 \
      --bfile ${scratch}/pgs/CT_disease_gwas_merged \
      --extract data/pgs/${phenotype}_disease_threshold_${threshold}.txt \
      --score data/gwas_results/${phenotype}_disease_combined.glm.logistic.hybrid 3 6 14 header no-mean-imputation \
      --memory 70000 \
      --out data/pgs/${phenotype}_disease_gwas_${threshold}_scores
  done
done
