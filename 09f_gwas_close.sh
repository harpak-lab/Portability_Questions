#!/bin/bash
#SBATCH -J GWAS_close
#SBATCH -o GWAS_close.o%j
#SBATCH -e GWAS_close.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

# Only run GWAS on the SNPs significant at 1e-5 in the original GWAS
for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC LDL Weight Triglycerides Cystatin_C Body_Fat_Perc
do
  $plink2 \
      --bfile ${scratch}/ukb_merged/ukb_merged \
      --keep data/ukb_populations/close_id.txt \
      --extract data/pgs/${phenotype}_threshold_1.txt \
      --pheno data/phenotypes/phenotypes.tsv \
      --pheno-name $phenotype \
      --require-pheno $phenotype \
      --covar data/ukb_merged/covar_close.tsv \
      --covar-name age sex_covar age_sq age_sex age_sq_sex $(printf "PC%i " $(seq 1 20)) \
      --vif 100000 \
      --memory 35000 \
      --glm hide-covar \
      --covar-variance-standardize \
      --no-input-missing-phenotype \
      --out data/gwas_results/${phenotype}_close
done
