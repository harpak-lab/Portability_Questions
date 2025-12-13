#!/bin/bash
#SBATCH -J extract_data_fields
#SBATCH -o extract_data_fields.o%j
#SBATCH -e extract_data_fields.o%j
#SBATCH -p normal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -A OTH21148
#SBATCH --mail-user=joyce.wang@utexas.edu
#SBATCH --mail-type=begin
#SBATCH --mail-type=end

file_handlers='/work/06568/joyce_w/stampede2/software/ukbconv'

# Sex chromosome aneuploidy
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22019 \
   -odata/extracted_data_fields/sex_chrom_aneuploidy

# Self-reported sex
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s31 \
   -odata/extracted_data_fields/reported_sex

# Sex determined from genotyping analysis
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22001 \
   -odata/extracted_data_fields/genetic_sex

# Outliers for heterozygosity or missing rate
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22027 \
   -odata/extracted_data_fields/outlier_heterozygosity

# Genotype missingness
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22005 \
   -odata/extracted_data_fields/genotype_missingness

# Used in PCA (for excluding close relatives)
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22020 \
   -odata/extracted_data_fields/used_in_pca

# Age
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s21022 \
   -odata/extracted_data_fields/age

# PC
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22009 \
   -odata/extracted_data_fields/pc

# Ethnic background
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s21000 \
   -odata/extracted_data_fields/ethnic_background

# Caucasian as determined by genetics
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22006 \
   -odata/extracted_data_fields/caucasian

# Genotype array used
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s22000 \
   -odata/extracted_data_fields/array_type

# Townsend deprivation index
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s189 \
   -odata/extracted_data_fields/townsend

# Total household income
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s738 \
   -odata/extracted_data_fields/income

# Educational attainment
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s6138 \
   -odata/extracted_data_fields/education

# ICD 10 main
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s41202 \
   -odata/extracted_data_fields/icd_10_main

# ICD 10 secondary
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s41204 \
   -odata/extracted_data_fields/icd_10_secondary
   
# Assessment center
$file_handlers/ukbconv \
   $file_handlers/ukb45020.enc_ukb \
   txt \
   -s54 \
   -odata/extracted_data_fields/assessment_center
