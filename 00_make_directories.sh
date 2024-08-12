#!/bin/bash
#SBATCH -J prepare_directories
#SBATCH -o prepare_directories.o%j
#SBATCH -e prepare_directories.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 00:05:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

mkdir data \
      data/extracted_data_fields \
      data/fst \
      data/genotype_calls \
      data/gwas_results \
      data/heterozygosity \
      data/info_filter \
      data/kgp_filtered \
      data/kgp_merged \
      data/kgp_meta \
      data/pca \
      data/phenotypes \
      data/pgs \
      data/pgs_pred \
      data/ukb_filtered \
      data/ukb_merged \
      data/ukb_populations \
      data/ukb_populations/bins \
      img

# For placing large data
scratch='/scratch/06568/joyce_w/'
mkdir ${scratch}/pgs_portability_questions \
      ${scratch}/pgs_portability_questions/data \
      ${scratch}/pgs_portability_questions/data/kgp \
      ${scratch}/pgs_portability_questions/data/kgp_filtered \
      ${scratch}/pgs_portability_questions/data/kgp_merged \
      ${scratch}/pgs_portability_questions/data/genotype_calls \
      ${scratch}/pgs_portability_questions/data/pgs \
      ${scratch}/pgs_portability_questions/data/ukb_filtered \
      ${scratch}/pgs_portability_questions/data/ukb_merged \
      ${scratch}/pgs_portability_questions/data/ukb_merged/geno_group
