#!/bin/bash
#SBATCH -J compare_genetic_variance
#SBATCH -o compare_genetic_variance.o%j
#SBATCH -e compare_genetic_variance.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 4:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'

# Directory
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'

# Calculate heterozygosity in the GWAS group
$plink \
  --bfile ${scratch}/ukb_merged/ukb_merged \
  --keep data/ukb_populations/wb_gwas_id.txt \
  --extract data/pgs/all_pgs_snps.txt \
  --freqx \
  --out data/heterozygosity/heterozygosity_wb

# Calculate heterozygosity in the "close" group
$plink \
  --bfile ${scratch}/ukb_merged/ukb_merged \
  --keep data/ukb_populations/close_id.txt \
  --extract data/pgs/all_pgs_snps.txt \
  --freqx \
  --out data/heterozygosity/heterozygosity_close

# Calculate heterozygosity in the "far" group
$plink \
  --bfile ${scratch}/ukb_merged/ukb_merged \
  --keep data/ukb_populations/far_id.txt \
  --extract data/pgs/all_pgs_snps.txt \
  --freqx \
  --out data/heterozygosity/heterozygosity_far

# Make the plots
Rscript 10b_heritability_plots.R

conda deactivate
