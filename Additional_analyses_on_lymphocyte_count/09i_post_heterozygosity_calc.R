library(dplyr)
library(tidyr)
library(readr)

# Group by magnitude
pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
          "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
          "Body_Fat_Perc")

# Get the SNPs used in PGS calculations
snp_lst <- list()

for(trait in pheno){
  snp <- read_tsv(paste0("data/pgs/", trait, "_threshold_1.txt"), 
                 col_names = F)
  colnames(snp) = "ID"
  
  beta <- read_tsv(paste0("data/gwas_results/", trait, "_combined.glm.linear"))
  snp <- snp %>% left_join(beta, by = "ID")
  snp$beta_sq <- (snp$BETA)^2
  
  snp <- snp[, c(1, 14)]
  
  snp_lst[[trait]] <- snp
}

# Stratify SNPs by into 3 equally sized bins
het <- c()
dir <- "data/heterozygosity/"

for(i in 1:250){
  print(i)
  temp <- read_tsv(paste0(dir, "heterozygosity_bin_", i, ".frqx"))
  for(trait in pheno){
    
    # Only look at SNPs used in PGS calculations
    temp2 <- temp %>% select(SNP, `C(HOM A1)`, `C(HET)`, `C(HOM A2)`)
    colnames(temp2)[1] <- "ID"
    temp2 <- snp_lst[[trait]] %>% left_join(temp2, by = "ID")
    
    temp2 <- temp2 %>% na.omit()
    temp2$beta_sq_size <- ntile(temp2$beta_sq, 3)
    
    # Calculate the heterozygosity of each SNP
    temp2$het <- temp2$`C(HET)` / ( temp2$`C(HOM A1)` +  temp2$`C(HET)` +  temp2$`C(HOM A2)`)
    # Take the mean of each stratum
    het_temp <- temp2 %>%
      dplyr::group_by(beta_sq_size) %>%
      dplyr::summarize(mean_het = mean(het))
    het <- c(het, het_temp$mean_het)
  }
}

heterozygosity <- cbind.data.frame(phenotype = rep(rep(pheno, each = 3), 250),
                                  weighted_pc_groups = rep(1:250, each = 15 * 3),
                                  size = rep(1:3, 250 * 15),
                                  mean_het = het)
heterozygosity <- heterozygosity %>%
  arrange(phenotype, weighted_pc_groups, size)

# Read in the file containing median PC distance of each bin
median_pc <- read_tsv('data/pgs_pred/group_pgs_df.tsv')
median_pc <- median_pc %>% distinct(weighted_pc_groups, median_pc, group_close_to_gwas)

heterozygosity <- heterozygosity %>% left_join(median_pc, by = "weighted_pc_groups")
heterozygosity %>% write_tsv("data/heterozygosity/mean_heterozygosity.tsv")
