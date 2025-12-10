#!/bin/bash
#SBATCH -J GWAS_array_center
#SBATCH -o GWAS_array_center.o%j
#SBATCH -e GWAS_array_center.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Creates covariate file
Rscript 14b_prepare_covariates_center_array.R

plink='/work/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC LDL Weight Triglycerides Cystatin_C Body_Fat_Perc
do
  for chromosome in $(seq 1 22);
  do
    $plink2 \
      --bfile ${scratch}/ukb_filtered/chr${chromosome}_filtered \
      --keep data/ukb_populations/wb_gwas_id.txt \
      --pheno data/phenotypes/phenotypes.tsv \
      --pheno-name $phenotype \
      --require-pheno $phenotype \
      --covar data/ukb_merged/covar_array_center.tsv \
      --covar-name age sex_covar age_sq age_sex age_sq_sex $(printf "PC%i " $(seq 1 20)) $(printf "UKBiLEVEAX_b%i " $(seq 11 -1 1)) $(printf "Batch_b%03i " $(seq 1 94)) Barts Birmingham Bristol Bury Cardiff Croydon Edinburgh Glasgow Hounslow Leeds Liverpool Manchester Middlesborough Newcastle Nottingham Oxford Reading Sheffield Stockport Stoke Swansea \
      --vif 100000 \
      --memory 70000 \
      --glm hide-covar \
      --covar-variance-standardize \
      --no-input-missing-phenotype \
      --out data/gwas_results/${phenotype}.chr${chromosome}_array_center
  done
done

conda deactivate
