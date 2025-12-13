#!/bin/bash
#SBATCH -J post-clump_regenie
#SBATCH -o post-clump_regenie.o%j
#SBATCH -e post-clump_regenie.o%j
#SBATCH -p normal
#SBATCH -N 3
#SBATCH -n 10
#SBATCH -t 48:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink="/work2/06568/joyce_w/stampede2/software/plink/plink/plink"

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

cat data/pgs/*_regenie_threshold_4.txt | uniq > data/pgs/all_pgs_regenie_snps.txt

for chromosome in $(seq 1 22);
do
  $plink \
    --bfile ${scratch}/ukb_filtered/chr${chromosome}_filtered \
    --extract data/pgs/all_pgs_regenie_snps.txt \
    --remove data/ukb_populations/wb_gwas_id.txt \
    --memory 70000 \
    --make-bed \
    --out ${scratch}/pgs/chr${chromosome}_regenie_temp

  printf "/scratch/06568/joyce_w/pgs_portability_questions/data/pgs/chr%s_regenie_temp\n" $chromosome >> data/pgs/pgs_regenie_merged_list.txt
done

# Combine extracted SNPs across chromosomes into a single Plink 1 file of all
# SNPs that meet the least extreme level of significance for any of the traits.
$plink \
  --merge-list data/pgs/pgs_regenie_merged_list.txt \
  --make-bed \
  --memory 70000 \
  --out ${scratch}/pgs/CT_regenie_merged

rm ${scratch}/pgs/chr*_temp*
