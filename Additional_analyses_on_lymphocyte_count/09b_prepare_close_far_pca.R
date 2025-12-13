library(tidyr)
library(dplyr)
library(readr)

# Read in the file containing genetic distance
pc_dist <- read_tsv("data/pgs_pred/group_non_pgs_df.tsv")
pc_dist <- pc_dist[, c(2, 12)]
pc_dist <- pc_dist %>% distinct(IID, pc_dist)

pc_dist <- cbind.data.frame(`#FID` = pc_dist$IID,
                            pc_dist)

# Get a list of IDs for individuals closer to the GWAS group
close <- pc_dist %>% filter(pc_dist <= 10) %>% select(-pc_dist)

# Get a list of IDs for individuals further away the GWAS group
far <- pc_dist %>% filter(pc_dist > 10) %>% select(-pc_dist)

close %>% write_tsv("data/ukb_populations/close_id.txt")
far %>% write_tsv("data/ukb_populations/far_id.txt")
