#!/bin/bash
#SBATCH -J ukb_kgp_pca
#SBATCH -o ukb_kgp_pca.o%j
#SBATCH -e ukb_kgp_pca.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 12:00:00
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

# Download the 1000 Genomes files from https://www.cog-genomics.org/plink/2.0/resources
wget -nd https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst
wget -nd https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz
wget -nd https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x

# Unzip the files
unzstd all_hg38.pgen.zst
mv 'all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz' all_hg38.pvar.zst
unzstd all_hg38.pvar.zst
mv 'hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x' all_hg38.psam
rm all_hg38*.zst

# Move the files
mv all_hg38* ${scratch}/kgp

# Download the meta file
wget -nd https://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/sample_info/integrated_call_samples.20130502.ALL.ped
mv integrated_call_samples.20130502.ALL.ped data/kgp_meta

# Download a list of SNPs used in UKB PCA
wget -nd biobank.ctsu.ox.ac.uk/ukb/ukb/auxdata/ukb_snp_qc.txt
mv ukb_snp_qc.txt data/pca

# Prepare files
Rscript 05i_prepare_files_for_ukb_kgp_pca.R

# Filter the 1000 Genomes file for CEU, CHB, and YRI individuals
$plink2 \
  --pfile ${scratch}/kgp/all_hg38 \
  --allow-extra-chr \
  --rm-dup 'exclude-mismatch' \
  --snps-only 'just-acgt' \
  --extract data/pca/snps_in_pca.txt \
  --max-alleles 2 \
  --keep data/kgp_meta/keep_id.tsv \
  --memory 70000 \
  --make-bed \
  --out ${scratch}/kgp_filtered/kgp_filtered

# Filter UKB file for SNPs used in UKB PCA
$plink \
    --bfile ${scratch}/ukb_merged/ukb_merged \
    --extract data/pca/snps_in_pca.txt \
    --make-bed \
    --out ${scratch}/ukb_merged/ukb_merged_snps_in_pca

# First merging attemp
$plink \
    --bfile ${scratch}/ukb_merged/ukb_merged_snps_in_pca \
    --bmerge ${scratch}/kgp_filtered/kgp_filtered \
    --allow-extra-chr \
    --merge-mode 5 \
    --make-bed \
    --out ${scratch}/kgp_filtered/ukb_kgp_merged

# Flip positions
$plink \
  --bfile ${scratch}/kgp_filtered/kgp_filtered \
  --allow-extra-chr \
  --flip ${scratch}/kgp_filtered/ukb_kgp_merged-merge.missnp \
  --make-bed \
  --out ${scratch}/kgp_merged/ukb_kgp_merged_source_2_trial

# Second merging attemp
$plink \
  --bfile ${scratch}/ukb_merged/ukb_merged_snps_in_pca \
  --allow-extra-chr \
  --bmerge ${scratch}/kgp_merged/ukb_kgp_merged_source_2_trial \
  --make-bed \
  --out ${scratch}/kgp_merged/ukb_kgp_merged_trial

# Remove problematic SNPs
$plink \
  --bfile ${scratch}/ukb_merged/ukb_merged_snps_in_pca \
  --exclude ${scratch}/kgp_merged/ukb_kgp_merged_trial-merge.missnp \
  --make-bed \
  --out ${scratch}/ukb_merged/merged_snps_in_pca_tmp

$plink \
  --bfile ${scratch}/kgp_merged/ukb_kgp_merged_source_2_trial \
  --allow-extra-chr \
  --exclude ${scratch}/kgp_merged/ukb_kgp_merged_trial-merge.missnp \
  --make-bed \
  --out ${scratch}/kgp_merged/ukb_kgp_merged_source_2_trial_tmp

# Final merging attemp
$plink \
  --bfile ${scratch}/ukb_merged/merged_snps_in_pca_tmp \
  --allow-extra-chr \
  --bmerge ${scratch}/kgp_merged/ukb_kgp_merged_source_2_trial_tmp \
  --make-bed \
  --out data/kgp_merged/ukb_kgp_merged

# Run PCA on the merged UKB and 1000 Genomes file
$plink2 \
  --memory 70000 \
  --bfile data/kgp_merged/ukb_kgp_merged \
  --pca approx 20 \
  --out data/pca/ukb_kgp_merged_pca

conda deactivate
