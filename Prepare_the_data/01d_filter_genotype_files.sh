#!/bin/bash
#SBATCH -J filter_genotype_files
#SBATCH -o filter_genotype_files.o%j
#SBATCH -e filter_genotype_files.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 4:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink="/work/06568/joyce_w/stampede2/software/plink/plink/plink"
plink2="/work/06568/joyce_w/stampede2/software/plink/plink2/plink2"

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

source /work/06568/joyce_w/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Remove the MHC region and extended region in strong LD with it
awk '$4>26477797 && $4<35448354' ${scratch}/genotype_calls/ukb_snp_chr6_v2.bim | awk '{print $2}' > data/mhc_variants.txt

for i in $(seq 1 22);
do

  # First write a list of SNPs to keep based on the GWAS set
  $plink2 \
    --bim ${scratch}/genotype_calls/ukb_snp_chr${i}_v2.bim \
    --bed ${scratch}/genotype_calls/ukb_cal_chr${i}_v2.bed \
    --fam ${scratch}/genotype_calls/ukb61666_cal_chr${i}_v2_s488180.fam \
    --keep data/ukb_populations/wb_gwas_id.txt \
    --exclude data/mhc_variants.txt \
    --hwe 1e-10 \
    --maf 0.0001 \
    --snps-only 'just-acgt' \
    --max-alleles 2 \
    --write-snplist \
    --out ${scratch}/ukb_filtered/chr${i}_to_keep

  $plink2 \
    --bim ${scratch}/genotype_calls/ukb_snp_chr${i}_v2.bim \
    --bed ${scratch}/genotype_calls/ukb_cal_chr${i}_v2.bed \
    --fam ${scratch}/genotype_calls/ukb61666_cal_chr${i}_v2_s488180.fam \
    --keep data/ukb_populations/keep_id.txt \
    --extract ${scratch}/ukb_filtered/chr${i}_to_keep.snplist \
    --make-bed \
    --out ${scratch}/ukb_filtered/chr${i}_filtered

  # Append the output file path to a new file (for merging them all below)
  printf "/scratch/06568/joyce_w/pgs_portability_questions/data/ukb_filtered/chr%s_filtered\n" $i >> data/ukb_merged/ukb_merged_list.txt
done

# Merge the files
$plink \
  --merge-list data/ukb_merged/ukb_merged_list.txt \
  --make-bed \
  --out ${scratch}/ukb_merged/ukb_merged
