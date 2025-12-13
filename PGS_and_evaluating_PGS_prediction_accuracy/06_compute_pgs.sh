#!/bin/bash
#SBATCH -J PGS
#SBATCH -o PGS.o%j
#SBATCH -e PGS.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 8:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e 

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC LDL Weight Triglycerides Cystatin_C Body_Fat_Perc
do
  for threshold in 0 1 2 3 4
  do
    # Score each individual using the SNPs below a p-value threshold and GWAS betas
    $plink2 \
      --bfile ${scratch}/pgs/CT_merged \
      --extract data/pgs/${phenotype}_threshold_${threshold}.txt \
      --score data/gwas_results/${phenotype}_combined.glm.linear 3 6 9 header no-mean-imputation \
      --memory 70000 \
      --out data/pgs/${phenotype}_${threshold}_scores
  done
done
