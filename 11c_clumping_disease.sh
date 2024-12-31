#!/bin/bash
#SBATCH -J CT_disease
#SBATCH -o CT_disease.o%j
#SBATCH -e CT_disease.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 10
#SBATCH -t 12:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e 

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
thresholds=(5e-8 1e-5 1e-4 1e-3 1e-2)

# Clump GWAS results 
for phenotype in Alzheimer T2D Asthma
do
  for chromosome in $(seq 1 22);
  do
    # Convert the GWAS output file from the Plink 2 to Plink 1 
    python 11d_convert_plink2_glm_to_plink1_disease.py \
      data/gwas_results/${phenotype}.chr${chromosome}_disease.${phenotype}.glm.logistic.hybrid \
      --output data/gwas_results/${phenotype}_disease.chr${chromosome}.${phenotype}.glm.assoc

    # Clump GWAS results using GWAS set
    $plink \
      --memory 70000 \
      --bfile ${scratch}/ukb_filtered/chr${chromosome}_filtered \
      --keep data/ukb_populations/wb_gwas_id.txt \
      --clump data/gwas_results/${phenotype}_disease.chr${chromosome}.${phenotype}.glm.assoc \
      --clump-p1 0.01 \
      --clump-r2 0.2 \
      --clump-kb 250 \
      --out data/gwas_results/${phenotype}_disease.chr${chromosome}.${phenotype}
  done

  # Combine clumped SNPs across chromosomes
  head -n 1 data/gwas_results/${phenotype}.chr1_disease.${phenotype}.clumped > data/gwas_results/${phenotype}_disease_combined.clumped 
  tail -n +2 -q data/gwas_results/${phenotype}.chr*_disease.${phenotype}.clumped >> data/gwas_results/${phenotype}_disease_combined.clumped 

  # Create files of SNPs meeting several p-value thresholds. Files numbered 0-4.
  for threshold in 0 1 2 3 4
  do
    # Further filter clumped SNPs using p-value thresholds (removes multiallelic SNPs)
    python 04c_filter_snps_for_pgs.py \
      data/gwas_results/${phenotype}_disease_combined.clumped \
      --threshold ${thresholds[$threshold]} \
      --output data/pgs/${phenotype}_disease_threshold_${threshold}.txt
  done

  # Create combined GWAS result files for each phenotype
  python 11e_combine_glm_threshold_4_disease.py \
    data/gwas_results/${phenotype}.chr*_disease.${phenotype}.glm.logistic.hybrid \
    --keep data/pgs/${phenotype}_disease_threshold_4.txt \
    --output data/gwas_results/${phenotype}_disease_combined.glm.logistic.hybrid
done

conda deactivate
