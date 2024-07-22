library(dplyr)
library(tidyr)
library(readr)

# Read in the file containing bins
bins <- read_tsv("data/pgs_pred/group_non_pgs_df.tsv")
bins <- bins[, c(1:2, 15)]
bins <- bins %>% distinct(`#FID`, IID, weighted_pc_groups)

# Write the IDs of each bin
for(i in 1:500){
  temp <- bins %>% filter(weighted_pc_groups == i)
  temp %>% 
    select(-weighted_pc_groups) %>% 
    write_tsv(paste0("data/ukb_populations/bins/bin_", i, "_id.txt"))
}