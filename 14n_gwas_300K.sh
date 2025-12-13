#!/bin/bash
#SBATCH -J GWAS_300K
#SBATCH -o GWAS_300K.o%j
#SBATCH -e GWAS_300K.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 10
#SBATCH -t 24:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink='/work/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

Rscript 14o_filter_individuals_300K.R

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC LDL Weight Triglycerides Cystatin_C Body_Fat_Perc
do
  for chromosome in $(seq 1 22);
  do
    $plink2 \
      --bfile ${scratch}/ukb_filtered/chr${chromosome}_filtered \
      --keep data/ukb_populations/wb_gwas_id_300K.txt \
      --pheno data/phenotypes/phenotypes.tsv \
      --pheno-name $phenotype \
      --require-pheno $phenotype \
      --covar data/ukb_merged/covar.tsv \
      --covar-name age sex age_sq age_sex age_sq_sex $(printf "PC%i " $(seq 1 20)) \
      --require-covar age sex age_sq age_sex age_sq_sex $(printf "PC%i " $(seq 1 20)) \
      --vif 100000 \
      --memory 70000 \
      --glm hide-covar \
      --covar-variance-standardize \
      --no-input-missing-phenotype \
      --out data/gwas_results/${phenotype}_300K.chr${chromosome}
  done
done

conda deactivate
