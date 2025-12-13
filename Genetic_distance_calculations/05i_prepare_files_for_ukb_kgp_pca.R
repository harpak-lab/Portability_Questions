library(readr)
library(tidyverse)

# Read in the 1000 Genomes meta file
meta <- read_tsv("data/kgp_meta/integrated_call_samples.20130502.ALL.ped")
# Filter for CEU, CHB, and YRI
meta <- meta %>% filter(Population %in% c("CEU", "CHB", "YRI"))
meta <- meta[, 1:2]
# Family IDs have to be set to 0 because the psam file doesn't contain these IDs
meta$`Family ID` = 0
# Write the IDs
meta %>% write_tsv("data/kgp_meta/keep_id.tsv")

# Get a list of SNPs used in UKB PCA
snp_pca <- read_table("data/pca/ukb_snp_qc.txt")
snp_pca <- snp_pca %>% filter(in_PCA == 1) %>%
  select(rs_id)
snp_pca %>% write.table("data/pca/snps_in_pca.txt", col.names = F, row.names = F, quote = F)
