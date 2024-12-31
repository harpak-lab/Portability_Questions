#!/bin/bash
#SBATCH -J GWAS_regenie
#SBATCH -o GWAS_regenie.o%j
#SBATCH -e GWAS_regenie.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate regenie_env

# For the phenotype and covariate files, the header of FID column starts with
# a "#", but regenie expect a header without a "#"
# Create a new version of both file without the "#"
sed '1s/^#//' data/phenotypes/phenotypes.tsv > data/phenotypes/phenotypes_regenie.tsv
sed '1s/^#//' data/ukb_merged/covar.tsv > data/ukb_merged/covar_regenie.tsv.tsv

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC LDL Weight Triglycerides Cystatin_C Body_Fat_Perc
do
  # Step 1
  regenie \
    --step 1 \
    --bed ${scratch}/ukb_filtered/chr${chromosome}_filtered \
    --keep data/ukb_populations/wb_gwas_id.txt \
    --phenoFile data/phenotypes/phenotypes_regenie.tsv \
    --phenoCol $phenotype \
    --covarFile data/ukb_merged/covar_regenie.tsv \
    --covarColList age,sex_covar,age_sq,age_sex,age_sq_sex,$(printf "PC%i," $(seq 1 20)) \
    --lowmem \
    --bsize 1000 \
    --threads 16 \
    --out data/gwas_results/${phenotype}_regenie_step_1
  
  # Step 2
  for chromosome in $(seq 1 22);
  do
    regenie \
      --step 2 \
      --bed ${scratch}/ukb_filtered/chr${chromosome}_filtered \
      --keep data/ukb_populations/wb_gwas_id.txt \
      --phenoFile data/phenotypes/phenotypes_regenie.tsv \
      --phenoCol $phenotype \
      --covarFile data/ukb_merged/covar_regenie.tsv \
      --covarColList age,sex_covar,age_sq,age_sex,age_sq_sex,$(printf "PC%i," $(seq 1 20)) \
      --lowmem \
      --bsize 500 \
      --threads 16 \
      --pred data/gwas_results/${phenotype}_regenie_step_1_pred.list \
      --out data/gwas_results/${phenotype}.chr${chromosome}_regenie
  done
done
