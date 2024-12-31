#!/bin/bash
#SBATCH -J PGS_prscs
#SBATCH -o PGS_prscs.o%j
#SBATCH -e PGS_prscs.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -t 16:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

set -e

plink='/work2/06568/joyce_w/stampede2/software/plink/plink/plink'
plink2='/work2/06568/joyce_w/stampede2/software/plink/plink2/plink2'
scratch='/scratch/06568/joyce_w/pgs_portability_questions/data'
prscs='/work2/06568/joyce_w/stampede2/software/PRScs'

source /work2/06568/joyce_w/stampede2/software/anaconda3/etc/profile.d/conda.sh
conda init bash
conda activate pgs

# Prepare summary statistics
Rscript 14m_prepare_sum_stats_prscs.R

conda deactivate

conda activate prscs

# Download and unzip the EUR LD reference panel
wget -P ${scratch}/ukb_filtered/ https://www.dropbox.com/s/mt6var0z96vb6fv/ldblk_1kg_eur.tar.gz
tar -zxvf ${scratch}/ukb_filtered/ldblk_1kg_eur.tar.gz -C ${scratch}/ukb_filtered/

for phenotype in BMI Lymphocyte Height Eosinophil MCH MCV Monocyte Platelet RBC WBC LDL Weight Triglycerides Cystatin_C Body_Fat_Perc
do
  for chromosome in $(seq 1 22);
  do
    python ${prscs}/PRScs.py \
      --ref_dir=${scratch}/ukb_filtered/ldblk_ukbb_eur \
      --bim_prefix=${scratch}/ukb_filtered/chr${chromosome}_filtered \
      --sst_file=data/gwas_results/${phenotype}.chr${chromosome}_prscs.${phenotype}.glm.linear \
      --n_gwas=350000 \
      --chrom=${chromosome} \
      --out_dir=data/pgs/${phenotype}_prscs_chr${chromosome}

    cat data/pgs/${phenotype}_prscs_chr${chromosome}_pst_eff_a1_b0.5_phiauto_chr${chromosome}.txt >> data/pgs/${phenotype}_prscs_combined.txt
  done
  
  $plink2 \
    --bfile ${scratch}/ukb_filtered/chr${chromosome}_filtered \
    --remove data/ukb_populations/wb_gwas_id.txt \
    --score data/pgs/${phenotype}_prscs_combined.txt 2 4 6 no-mean-imputation \
    --memory 10000 \
    --threads 16 \
    --out data/pgs/${phenotype}_prscs_scores
done

conda deactivate
