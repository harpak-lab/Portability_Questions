#!/bin/bash
#SBATCH -J prepare_for_ma_counts
#SBATCH -o prepare_for_ma_counts.o%j
#SBATCH -e prepare_for_ma_counts.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 24:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

# Export the genotypes
$plink2 \
  --bfile ${scratch}/ukb_merged/ukb_merged \
  --extract data/pgs/all_pgs_snps.txt
  --export A \
  --out ${scratch}/ukb_merged/ukb_merged_geno

# Subset the genotype file so it's easier to process
group=1
low=$((($group-1)*1000+7))
high=$(($group*1000+6))

while [ "$group" -le 169 ]
do
  echo "$group, $low, $high"
  cut -f1-2,${low}-${high} ${scratch}/ukb_merged/ukb_merged_geno.raw > ${scratch}/ukb_merged/geno_group/group_$group.txt

  low=$(($low+1000))
  high=$(($high+1000))
  group=$(($group+1))
done

echo "$group, $low, $high"

input_file="/scratch/06568/joyce_w/pgs_portability_questions/data/ukb_merged/ukb_merged_geno.raw"
max=$(awk '{print NF; exit}' "$input_file")

cut -f1-2,${low}-${max} ${scratch}/ukb_merged/ukb_merged_geno.raw > ${scratch}/ukb_merged/geno_group/group_$group.txt
