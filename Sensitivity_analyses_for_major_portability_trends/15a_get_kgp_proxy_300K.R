library(tidyr)
library(dplyr)
library(readr)

# Read in PCA and meta files
pc <- read_tsv("data/pca/ukb_kgp_merged_pca.eigenvec")
pc <- pc[, 2:22]
meta <- read_tsv("data/kgp_meta/integrated_call_samples.20130502.ALL.ped")
meta <- meta %>% select(`Individual ID`, Population)
colnames(meta) <- c("IID", "Pop")
pc <- left_join(pc, meta, by = "IID")
pc <- cbind.data.frame(Pop = pc$Pop, pc[, 1:21])

# Get the centroid of each of the population
ceu_pc <- pc %>% filter(Pop == "CEU")
ceu_centroid <- colMeans(ceu_pc[,3:22])

chb_pc <- pc %>% filter(Pop == "CHB")
chb_centroid <- colMeans(chb_pc[,3:22])

yri_pc <- pc %>% filter(Pop == "YRI")
yri_centroid <- colMeans(yri_pc[,3:22])

# Read in the eigenvalues
eigenval <- read.table("data/pca/ukb_kgp_merged_pca.eigenval")
eigenval <- eigenval$V1

# Calculate weighted Euclidean distance of each individual to each of the centroids
get_pc_dist <- function(pc_val, centroid, var_explained, last_pc_dist_sq = 0){
  dist_sq <- (pc_val - centroid) ^ 2 * var_explained + last_pc_dist_sq
  return(dist_sq)
}

pc_dist_ceu_sq <- pc[, c(1,2)]
for(i in 1:20){
  if(i == 1){
    pc_dist_ceu_sq[["PC1_dist"]] <- get_pc_dist(pc[["PC1"]], ceu_centroid[["PC1"]],
                                               eigenval[1])
  } else{
    pc_dist_ceu_sq[[paste0("PC", i, "_dist")]] <- get_pc_dist(pc[[paste0("PC", i)]], 
                                                             ceu_centroid[[paste0("PC", i)]],
                                                             eigenval[i],
                                                             pc_dist_ceu_sq[[paste0("PC", i - 1, "_dist")]])
  }
}

pc_dist_ceu <- as.data.frame(apply(pc_dist_ceu_sq[, c(3:22)], 2, sqrt))
pc_dist_ceu <- cbind.data.frame(pc[, c(1:2)], pc_dist_ceu)

pc_dist_chb_sq <- pc[, c(1,2)]
for(i in 1:20){
  if(i == 1){
    pc_dist_chb_sq[["PC1_dist"]] <- get_pc_dist(pc[["PC1"]], chb_centroid[["PC1"]],
                                               eigenval[1])
  } else{
    pc_dist_chb_sq[[paste0("PC", i, "_dist")]] <- get_pc_dist(pc[[paste0("PC", i)]], 
                                                             chb_centroid[[paste0("PC", i)]],
                                                             eigenval[i],
                                                             pc_dist_chb_sq[[paste0("PC", i - 1, "_dist")]])
  }
}

pc_dist_chb <- as.data.frame(apply(pc_dist_chb_sq[, c(3:22)], 2, sqrt))
pc_dist_chb <- cbind.data.frame(pc[, c(1:2)], pc_dist_chb)

pc_dist_yri_sq <- pc[, c(1,2)]
for(i in 1:20){
  if(i == 1){
    pc_dist_yri_sq[["PC1_dist"]] <- get_pc_dist(pc[["PC1"]], yri_centroid[["PC1"]],
                                               eigenval[1])
  } else{
    pc_dist_yri_sq[[paste0("PC", i, "_dist")]] <- get_pc_dist(pc[[paste0("PC", i)]], 
                                                             yri_centroid[[paste0("PC", i)]],
                                                             eigenval[i],
                                                             pc_dist_yri_sq[[paste0("PC", i - 1, "_dist")]])
  }
}

pc_dist_yri <- as.data.frame(apply(pc_dist_yri_sq[, c(3:22)], 2, sqrt))
pc_dist_yri <- cbind.data.frame(pc[, c(1:2)], pc_dist_yri)


# Identify the prediction individual closest to each of the centroids
pc_dist_ceu_closest <- pc_dist_ceu %>% 
  filter(is.na(Pop)) %>%
  filter(PC20_dist == min(PC20_dist, na.rm = T))

pc_dist_chb_closest <- pc_dist_chb %>% 
  filter(is.na(Pop)) %>%
  filter(PC20_dist == min(PC20_dist, na.rm = T))

pc_dist_yri_closest <- pc_dist_yri %>% 
  filter(is.na(Pop)) %>%
  filter(PC20_dist == min(PC20_dist, na.rm = T))

# Use the genetic distance of the individual closest to the 3 1000 Genomes centroids as a proxy
best_pc_dist_gwas_std <- read_tsv("data/pca/pc_dist_best_gwas_300K_std.tsv")
best_pc_dist_pred_std <- read_tsv("data/pca/pc_dist_best_pred_300K_std.tsv")
best_pc_dist <- cbind.data.frame(best_pc_dist_gwas_std, best_pc_dist_pred_std)

proxy <- c(best_pc_dist$pc_dist[best_pc_dist$IID == pc_dist_ceu_closest$IID],
          best_pc_dist$pc_dist[best_pc_dist$IID == pc_dist_chb_closest$IID],
          best_pc_dist$pc_dist[best_pc_dist$IID == pc_dist_yri_closest$IID])

proxy <- cbind.data.frame(Pop = c("CEU", "CHB", "YRI"), pc_dist = proxy)

proxy %>% write_tsv("data/kgp_merged/pca_proxy_300K.tsv")
