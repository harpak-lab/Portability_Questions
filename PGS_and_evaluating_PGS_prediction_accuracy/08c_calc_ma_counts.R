library(tidyr)
library(dplyr)
library(readr)

pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
          "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
          "Body_Fat_Perc")

# Read in the GWAS IDs
gwas_id <- read.table("data/ukb_populations/wb_gwas_id.txt",
                     col.names = c("FID", "IID"))

`%notin%` = Negate(`%in%`)

# Get the SNPs used in PGS calculations
snp_lst <- list()

for(trait in pheno){
  snp <- read_tsv(paste0("data/pgs/", trait, "_threshold_1.txt"), 
                 col_names = F)
  colnames(snp) <- "ID"
  
  beta <- read_tsv(paste0("data/gwas_results/", trait, "_combined.glm.linear"))
  snp <- snp %>% left_join(beta, by = "ID")
  snp$beta_sq = (snp$BETA)^2
  
  snp <- snp[, c(1, 14)]
  
  snp_lst[[trait]] <- snp
}

# Read in genotypes
geno <- read_tsv("/scratch/06568/joyce_w/pgs_portability_questions/data/ukb_merged/geno_group/group_1.txt")
geno_gwas <- geno %>% filter(IID %in% gwas_id$IID)

# Count the # of minor alleles with respect to the GWAS group
count_0 <- colSums(geno_gwas[, 3:ncol(geno)] == 0, na.rm = T)
count_2 <- colSums(geno_gwas[, 3:ncol(geno)] == 2, na.rm = T)

count_df <- cbind.data.frame(id = colnames(geno_gwas)[3:ncol(geno)], 
                             count_0 = count_0, 
                             count_2 = count_2)

count_df$minor <- ifelse(count_df$count_0 < count_df$count_2, 0, 2)

geno_test <- geno %>% filter(IID %notin% gwas_id$IID)

minor_count <- geno_test[,2]

for(i in 3:ncol(geno)){
  minor_count[[colnames(geno_test)[i]]] <- ifelse(geno_test[[i]] == 1, 1,
                                                  ifelse(geno_test[[i]] == count_df$minor[i - 2], 2, 0))
}

minor_count_strat <- list()

temp <- minor_count %>% 
  pivot_longer(cols = 2:ncol(minor_count), names_to = "ID", values_to = "count")
temp$ID <- sub("_[^_]+$", "", temp$ID)

for(trait in pheno){
  minor_count_strat_df <- minor_count[, 1]
  minor_count_strat_df$low <- 0
  minor_count_strat_df$medium <- 0
  minor_count_strat_df$high <- 0
  
  temp2 <- temp %>% filter(ID %in% snp_lst[[trait]]$ID)
  temp2 <- temp2 %>% left_join(snp_lst[[trait]], by = "ID")
  
  # Separate the SNPs into 3 tiless based on the sizes of beta^2
  tile <- cbind.data.frame(ID = snp_lst[[trait]]$ID,
                           tile = ntile(snp_lst[[trait]]$beta_sq, 3))
  
  snps_low <- temp2 %>% filter(ID %in% tile$ID[tile$tile == 1])
  snps_medium <- temp2 %>% filter(ID %in% tile$ID[tile$tile == 2])
  snps_high <- temp2 %>% filter(ID %in% tile$ID[tile$tile == 3])
  
  if(nrow(snps_low) != 0){
    snps_low <- snps_low %>% 
      group_by(IID) %>%
      summarize(count_sum = sum(count, na.rm = T))
    
    minor_count_strat_df <- minor_count_strat_df %>% 
      left_join(snps_low, by = "IID")
    minor_count_strat_df$low <- minor_count_strat_df$low + minor_count_strat_df$count_sum
    minor_count_strat_df <- minor_count_strat_df %>% select(-count_sum)
  }
  
  if(nrow(snps_medium) != 0){
    snps_medium <- snps_medium %>% 
      group_by(IID) %>%
      summarize(count_sum = sum(count, na.rm = T))
    
    minor_count_strat_df <- minor_count_strat_df %>% 
      left_join(snps_medium, by = "IID")
    minor_count_strat_df$medium <- minor_count_strat_df$medium + minor_count_strat_df$count_sum
    minor_count_strat_df <- minor_count_strat_df %>% select(-count_sum)
  }
  
  if(nrow(snps_high) != 0){
    snps_high <- snps_high %>% 
      group_by(IID) %>%
      summarize(count_sum = sum(count, na.rm = T))
    
    minor_count_strat_df <- minor_count_strat_df %>% 
      left_join(snps_high, by = "IID")
    minor_count_strat_df$high <- minor_count_strat_df$high + minor_count_strat_df$count_sum
    minor_count_strat_df <- minor_count_strat_df %>% select(-count_sum)
  }
  
  minor_count_strat[[trait]] <- minor_count_strat_df
}

# Repeat for the other groups
for(group in 2:170){
  print(group)
  
  # Read in genotypes
  geno <- read_tsv(paste0("/scratch/06568/joyce_w/pgs_portability_questions/data/ukb_merged/geno_group/group_", group, ".txt"))
  geno_gwas <- geno %>% filter(IID %in% gwas_id$IID)
  
  # Count the # of minor alleles with respect to the GWAS group
  count_0 <- colSums(geno_gwas[, 3:ncol(geno)] == 0, na.rm = T)
  count_2 <- colSums(geno_gwas[, 3:ncol(geno)] == 2, na.rm = T)
  
  count_df <- cbind.data.frame(id = colnames(geno_gwas)[3:ncol(geno)], 
                               count_0 = count_0, 
                               count_2 = count_2)
  
  count_df$minor <- ifelse(count_df$count_0 < count_df$count_2, 0, 2)
  
  geno_test <- geno %>% filter(IID %notin% gwas_id$IID)
  
  minor_count <- geno_test[,2]
  
  for(i in 3:ncol(geno)){
    minor_count[[colnames(geno_test)[i]]] <- ifelse(geno_test[[i]] == 1, 1,
                                                    ifelse(geno_test[[i]] == count_df$minor[i - 2], 2, 0))
  }
  
  temp <- minor_count %>% pivot_longer(cols = 2:ncol(minor_count), names_to = "ID", values_to = "count")
  temp$ID <- sub("_[^_]+$", "", temp$ID)
  
  for(trait in pheno){
    minor_count_strat_df <- minor_count_strat[[trait]]
    
    temp2 <- temp %>% filter(ID %in% snp_lst[[trait]]$ID)
    temp2 <- temp2 %>% left_join(snp_lst[[trait]], by = "ID")
    
    # Separate the SNPs into 3 tiless based on the sizes of beta^2
    tile <- cbind.data.frame(ID = snp_lst[[trait]]$ID,
                            tile = ntile(snp_lst[[trait]]$beta_sq, 3))
    
    snps_low <- temp2 %>% filter(ID %in% tile$ID[tile$tile == 1])
    snps_medium <- temp2 %>% filter(ID %in% tile$ID[tile$tile == 2])
    snps_high <- temp2 %>% filter(ID %in% tile$ID[tile$tile == 3])
    
    if(nrow(snps_low) != 0){
      snps_low <- snps_low %>% 
        group_by(IID) %>%
        summarize(count_sum = sum(count, na.rm = T))
      
      minor_count_strat_df <- minor_count_strat_df %>% 
        left_join(snps_low, by = "IID")
      minor_count_strat_df$low <- minor_count_strat_df$low + minor_count_strat_df$count_sum
      minor_count_strat_df <- minor_count_strat_df %>% select(-count_sum)
    }
    
    if(nrow(snps_medium) != 0){
      snps_medium <- snps_medium %>% 
        group_by(IID) %>%
        summarize(count_sum = sum(count, na.rm = T))
      
      minor_count_strat_df <- minor_count_strat_df %>% 
        left_join(snps_medium, by = "IID")
      minor_count_strat_df$medium <- minor_count_strat_df$medium + minor_count_strat_df$count_sum
      minor_count_strat_df <- minor_count_strat_df %>% select(-count_sum)
    }
    
    if(nrow(snps_high) != 0){
      snps_high <- snps_high %>% 
        group_by(IID) %>%
        summarize(count_sum = sum(count, na.rm = T))
      
      minor_count_strat_df <- minor_count_strat_df %>% 
        left_join(snps_high, by = "IID")
      minor_count_strat_df$high <- minor_count_strat_df$high + minor_count_strat_df$count_sum
      minor_count_strat_df <- minor_count_strat_df %>% select(-count_sum)
    }
    
    minor_count_strat[[trait]] <- minor_count_strat_df
  }
}

# Convert the list into a data frame
minor_allele_count <- minor_count_strat[[pheno[1]]]
minor_allele_count <- cbind.data.frame(phenotype = pheno[1], 
                                       minor_allele_count)

for(trait in pheno[2:15]){
  minor_allele_count_temp <- minor_count_strat[[trait]]
  minor_allele_count_temp <- cbind.data.frame(phenotype = trait, 
                                              minor_allele_count_temp)
  minor_allele_count <- rbind.data.frame(minor_allele_count, minor_allele_count_temp)
}

minor_allele_count %>% write_tsv("data/ukb_merged/minor_allele_count.tsv")
