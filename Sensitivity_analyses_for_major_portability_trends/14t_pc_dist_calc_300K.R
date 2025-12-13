library(tidyr)
library(dplyr)
library(readr)
library(data.table)

# Get IDs
keep_id <- read.table("data/ukb_populations/keep_id.txt")
colnames(keep_id) = c("#FID", "IID")

gwas_id <- read.table("data/ukb_populations/wb_gwas_id_300K.txt")
colnames(gwas_id) = c("#FID", "IID")

# Read PCs from the file containing covariates
pc <- read_tsv('data/ukb_merged/covar.tsv')
pc <- pc[, -c(1, 3:7, 48)]
pc <- left_join(keep_id, pc, by = "IID")

# Get variance explained by each PC
pc_var <- apply(pc[, c(3:42)], 2, var)
# Standardize by the sum of variance explained
pc_var_std <- pc_var / sum(pc_var)

# Get the coordinates of the GWAS centroid
gwas_pc <- pc %>% filter(IID %in% gwas_id$IID)
gwas_centroid <- colMeans(gwas_pc[, 3:42])

# Calculate the squared Euclidean distance, weighted by variance explained
# Each value is the distance between current PC and the GWAS centroid, squared, multiplied by
# the variance explained, and with the value from the previous PC added
get_pc_dist = function(pc_val, gwas_centroid, var_explained, last_pc_dist_sq = 0){
  dist_sq = (pc_val - gwas_centroid) ^ 2 * var_explained + last_pc_dist_sq
  return(dist_sq)
}

pc_dist_sq = pc[, c(1,2)]
for(i in 1:40){
  if(i == 1){
    pc_dist_sq[["PC1_dist"]] <- get_pc_dist(pc[["PC1"]], gwas_centroid[["PC1"]], 
                                           pc_var_std[["PC1"]])
  } else{
    pc_dist_sq[[paste0("PC", i, "_dist")]] <- get_pc_dist(pc[[paste0("PC", i)]], 
                                                         gwas_centroid[[paste0("PC", i)]],
                                                         pc_var_std[[paste0("PC", i)]],
                                                         pc_dist_sq[[paste0("PC", i - 1, "_dist")]])
  }
}

# Take the sq root of all the values
pc_dist <- as.data.frame(apply(pc_dist_sq[, c(3:42)], 2, sqrt))
pc_dist <- cbind.data.frame(pc[, c(1:2)], pc_dist)
pc_dist %>% write_tsv("data/pca/pc_dist_300K.tsv")

# Get the 95th percentile of the PC distance of the GWAS individuals based on the first 10 PCs
gwas_pc_dist <- pc_dist %>% filter(pc_dist$IID %in% gwas_id$IID)
perc_95 <- quantile(gwas_pc_dist$PC10_dist, 0.95) 

# Prediction samples who are further away from the GWAS set
# Must have a PC distance > 95th percentile of the PC distance of the GWAS individuals based on the first 10 PCs
pc_dist_further <- pc_dist %>% 
  filter(pc_dist$PC10_dist > perc_95[[1]] & pc_dist$IID %notin% gwas_id$IID)

# Sample 10000 individuals to calculate Fst for
set.seed(3)
pc_dist_further_rand <- sample(pc_dist_further$IID, 10000)
pc_dist_further_rand <- sort(pc_dist_further_rand)
pc_dist_further_rand %>% 
  as.data.frame() %>% 
  write_tsv("data/ukb_populations/pc_futher_fst_id_300K.txt",
            col_names = F)
