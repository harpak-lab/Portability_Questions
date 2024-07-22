library(tidyr)
library(dplyr)
library(readr)
library(data.table)

# Read Fst and PC distance files
fst <- read_tsv('data/fst/final_fst.tsv')
pc_dist <- read_tsv("data/pca/pc_dist.tsv")
pc_dist_further_rand_id <- read_tsv("data/ukb_populations/pc_futher_fst_id.txt",
                                   col_names = "IID")
pc_dist_further_rand_id <- pc_dist_further_rand_id$IID
pc_dist_further_rand <- pc_dist %>% filter(pc_dist$IID %in% pc_dist_further_rand_id)
pc_dist_further_rand <- pc_dist_further_rand %>% left_join(fst, by = "IID")

# Correlate Fst with the PC distance calculated using different numbers of PCs
corr = c()
for(i in 1:40){
  correlation <- cor(pc_dist_further_rand[[paste0("PC", i, "_dist")]], 
                    pc_dist_further_rand$Weighted_Fst, use = "complete.obs")
  corr <- c(corr, correlation)
}

# Find the number of PCs that produces the highest correlation with Fst
max_corr <- max(corr)
max_corr_num = 0
for (i in 1:40){
  if(corr[i] == max_corr){
    max_corr_num <- i
  }
}

# Use that number of PCs
best_pc_dist <- cbind.data.frame(pc_dist[, c(1:2)], 
                                pc_dist = pc_dist[[paste0("PC", max_corr_num, "_dist")]])

best_pc_dist %>% write_tsv("data/pca/pc_dist_best.tsv")

# Separete the GWAS and prediction set
wb_gwas <- read_delim('data/ukb_populations/wb_gwas_id.txt', delim = ' ', trim_ws = T)

`%notin%` = Negate(`%in%`)

best_pc_dist_gwas <- best_pc_dist %>% filter(best_pc_dist$IID %in% wb_gwas$IID)
best_pc_dist_pred <- best_pc_dist %>% filter(best_pc_dist$IID %notin% wb_gwas$IID)

best_pc_dist_gwas %>% write_tsv("data/pca/pc_dist_best_gwas.tsv")
best_pc_dist_pred %>% write_tsv("data/pca/pc_dist_best_pred.tsv")

# Standardize the PC distance according to the mean PC distance of the GWAS set
best_pc_dist_gwas_std = cbind.data.frame(best_pc_dist_gwas[, c(1:2)],
                                         pc_dist = best_pc_dist_gwas$pc_dist / mean(best_pc_dist_gwas$pc_dist))

best_pc_dist_pred_std = cbind.data.frame(best_pc_dist_pred[, c(1:2)],
                                         pc_dist = best_pc_dist_pred$pc_dist / mean(best_pc_dist_gwas$pc_dist))

best_pc_dist_gwas_std %>% write_tsv("data/pca/pc_dist_best_gwas_std.tsv")
best_pc_dist_pred_std %>% write_tsv("data/pca/pc_dist_best_pred_std.tsv")
