library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(stringr)

# Phenotypes
pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
          "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
          "Body_Fat_Perc")

# Read in PGS for threshold 1 (p < 1e-5)
pgs_df <- data.frame()
for(trait in pheno){
  temp <- read_tsv(paste0("data/pgs/", trait, "_1_scores.sscore"))
  temp$pgs <- temp$SCORE1_AVG
  temp <- temp %>% select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG))
  temp$phenotype <- rep(trait, nrow(temp))
  pgs_df <- rbind.data.frame(pgs_df, temp)
  rm(temp)
}

# Read in group level file for covariates
group_non_pgs_df <- read_tsv("data/pgs_pred/group_non_pgs_df.tsv")

pgs_df <- left_join(group_non_pgs_df, pgs_df, by = c("#FID", "IID", "phenotype"))
pgs_df <- pgs_df %>% select(`#FID`, IID, sex, age, age_sq, age_sex, age_sq_sex, array_type,
                           pc_dist, pgs, phenotype, phenotype_value, weighted_pc_groups, group_close_to_gwas)

pgs_df <- pgs_df %>% mutate(
  sex_pc_dist = sex * pc_dist,
  pc_dist_sq = pc_dist ^ 2,
  pc_dist_3 = pc_dist ^ 3,
  pc_dist_4 = pc_dist ^ 4,
  pc_dist_5 = pc_dist ^ 5,
  pc_dist_6 = pc_dist ^ 6,
  pc_dist_7 = pc_dist ^ 7,
  pc_dist_8 = pc_dist ^ 8,
  pc_dist_9 = pc_dist ^ 9,
  pc_dist_10 = pc_dist ^ 10,
  pc_dist_11 = pc_dist ^ 11,
  pc_dist_12 = pc_dist ^ 12,
  pc_dist_13 = pc_dist ^ 13,
  pc_dist_14 = pc_dist ^ 14,
  pc_dist_15 = pc_dist ^ 15,
  pc_dist_16 = pc_dist ^ 16,
  pc_dist_17 = pc_dist ^ 17,
  pc_dist_18 = pc_dist ^ 18,
  pc_dist_19 = pc_dist ^ 19,
  pc_dist_20 = pc_dist ^ 20
)

# Fit a linear model of Y ~ covariates in each bin independently, 
# startng with the first bin
temp <- pgs_df %>% filter(phenotype == pheno[1] & weighted_pc_groups == 1)
# Only use array_type if both array types have more than 1 sample in this bin
# (or lm doesn't run)
if(length(temp$array_type[temp$array_type == "Axiom"]) > 1 &
   length(temp$array_type[temp$array_type == "BiLEVE"]) > 1){
  Y_mod <- temp %>%
    lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex + array_type, data = .)
  
  # Apply the model onto this weighted_pc_groups
  temp$Y_hat <- unname(predict(Y_mod, newdata = temp))
  
} else{
  Y_mod <- temp %>%
    lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex, data = .)
  
  # Apply the model onto this weighted_pc_groups
  temp$Y_hat <- unname(predict(Y_mod, newdata = temp))
}

# Repeat for the other bins
for(i in 2:500){
  temp2 <- pgs_df %>% filter(phenotype == pheno[1] & weighted_pc_groups == i)
  if(length(temp2$array_type[temp2$array_type == "Axiom"]) > 1 &
     length(temp2$array_type[temp2$array_type == "BiLEVE"]) > 1){
    Y_mod <- temp2 %>%
      lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex + array_type, data = .)
    
    temp2$Y_hat = unname(predict(Y_mod, newdata = temp2))
    
  } else{
    Y_mod <- temp2 %>%
      lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex, data = .)
    
    temp2$Y_hat <- unname(predict(Y_mod, newdata = temp2))
  }
  temp <- rbind.data.frame(temp, temp2)
}

# Do the same for the other phenotypes
for(trait in pheno[2:15]){
  print(trait)
  for(i in 1:500){
    temp2 <- pgs_df %>% filter(phenotype == trait & weighted_pc_groups == i)
    if(length(temp2$array_type[temp2$array_type == "Axiom"]) > 1 &
       length(temp2$array_type[temp2$array_type == "BiLEVE"]) > 1){
      Y_mod <- temp2 %>%
        lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex + array_type, data = .)
      
      temp2$Y_hat = unname(predict(Y_mod, newdata = temp2))
      
    } else{
      Y_mod <- temp2 %>%
        lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex, data = .)
      
      temp2$Y_hat <- unname(predict(Y_mod, newdata = temp2))
    }
    temp <- rbind.data.frame(temp, temp2)
  }
}

pgs_df <- temp %>% arrange(IID)

# X = Y - Y_hat
pgs_df$X <- pgs_df$phenotype_value - pgs_df$Y_hat

# Do X ~ gen_dist_polynomial + sex + sex_gen_dist across all bins
X_mod <- pgs_df %>%
  group_by(phenotype) %>%
  do(lm_X = lm(X ~ pc_dist + pc_dist_sq + pc_dist_3 + pc_dist_4 + pc_dist_5 + pc_dist_6 +
                 pc_dist_7 + pc_dist_8 + pc_dist_9 + pc_dist_10 + pc_dist_11 + pc_dist_12 +
                 pc_dist_13 + pc_dist_14 + pc_dist_15 + pc_dist_16 + pc_dist_17 + pc_dist_18 +
                 pc_dist_19 + pc_dist_20 + sex + sex_pc_dist, data = .))

# Apply this model
temp <- pgs_df %>% filter(phenotype == pheno[1])
temp$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df %>% filter(phenotype == trait)
  temp2$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df <- temp %>% arrange(IID)

# Z = X - X_hat
pgs_df$Z <- pgs_df$X - pgs_df$X_hat

# Z ~ PGS
Z_mod <- pgs_df %>%
  group_by(phenotype) %>%
  do(lm_Z = lm(Z ~ pgs, data = .))

# Apply this model
temp <- pgs_df %>% filter(phenotype == pheno[1])
temp$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df %>% filter(phenotype == trait)
  temp2$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df <- temp %>% arrange(IID)

# Prediction error = (Z - Z_hat) ^ 2
pgs_df$pred_error = (pgs_df$Z - pgs_df$Z_hat) ^ 2

pgs_df <- pgs_df %>% arrange(IID, phenotype)

pgs_df %>% write_tsv("data/pgs_pred/ind_pgs_df.tsv")