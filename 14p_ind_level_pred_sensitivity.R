library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(stringr)

# Phenotypes
pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
           "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
           "Body_Fat_Perc")

# Read in group level file for covariates
group_non_pgs_df <- read_tsv("data/pgs_pred/group_non_pgs_df.tsv")

# GWAS with assessment center and genotype array
# Read in PGS for threshold 1 (p < 1e-5)
pgs_df_array_center <- data.frame()
for(trait in pheno){
  temp <- read_tsv(paste0("data/pgs/", trait, "_array_center_1_scores.sscore"))
  temp$pgs <- temp$SCORE1_AVG
  temp <- temp %>% select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG))
  temp$phenotype <- rep(trait, nrow(temp))
  pgs_df_array_center <- rbind.data.frame(pgs_df_array_center, temp)
  rm(temp)
}

pgs_df_array_center <- left_join(group_non_pgs_df, pgs_df_array_center, by = c("#FID", "IID", "phenotype"))
pgs_df_array_center <- pgs_df_array_center %>% select(`#FID`, IID, sex, age, age_sq, age_sex, age_sq_sex, array_type,
                            pc_dist, pgs, phenotype, phenotype_value, weighted_pc_groups, group_close_to_gwas)

pgs_df_array_center <- pgs_df_array_center %>% mutate(
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
temp <- pgs_df_array_center %>% filter(phenotype == pheno[1] & weighted_pc_groups == 1)
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
for(i in 2:250){
  temp2 <- pgs_df_array_center %>% filter(phenotype == pheno[1] & weighted_pc_groups == i)
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
  for(i in 1:250){
    temp2 <- pgs_df_array_center %>% filter(phenotype == trait & weighted_pc_groups == i)
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

pgs_df_array_center <- temp %>% arrange(IID)

# X = Y - Y_hat
pgs_df_array_center$X <- pgs_df_array_center$phenotype_value - pgs_df_array_center$Y_hat

# Do X ~ gen_dist_polynomial + sex + sex_gen_dist across all bins
X_mod <- pgs_df_array_center %>%
  group_by(phenotype) %>%
  do(lm_X = lm(X ~ pc_dist + pc_dist_sq + pc_dist_3 + pc_dist_4 + pc_dist_5 + pc_dist_6 +
                 pc_dist_7 + pc_dist_8 + pc_dist_9 + pc_dist_10 + pc_dist_11 + pc_dist_12 +
                 pc_dist_13 + pc_dist_14 + pc_dist_15 + pc_dist_16 + pc_dist_17 + pc_dist_18 +
                 pc_dist_19 + pc_dist_20 + sex + sex_pc_dist, data = .))

# Apply this model
temp <- pgs_df_array_center %>% filter(phenotype == pheno[1])
temp$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_array_center %>% filter(phenotype == trait)
  temp2$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_array_center <- temp %>% arrange(IID)

# Z = X - X_hat
pgs_df_array_center$Z <- pgs_df_array_center$X - pgs_df_array_center$X_hat

# Z ~ PGS
Z_mod <- pgs_df_array_center %>%
  group_by(phenotype) %>%
  do(lm_Z = lm(Z ~ pgs, data = .))

# Apply this model
temp <- pgs_df_array_center %>% filter(phenotype == pheno[1])
temp$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_array_center %>% filter(phenotype == trait)
  temp2$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_array_center <- temp %>% arrange(IID)

# Prediction error = (Z - Z_hat) ^ 2
pgs_df_array_center$pred_error = (pgs_df_array_center$Z - pgs_df_array_center$Z_hat) ^ 2

pgs_df_array_center <- pgs_df_array_center %>% arrange(IID, phenotype)

pgs_df_array_center %>% write_tsv("data/pgs_pred/ind_pgs_df_array_center.tsv")


# GWAS with regenie
# Read in PGS for threshold 1 (p < 1e-5)
pgs_df_regenie <- data.frame()
for(trait in pheno){
  temp <- read_tsv(paste0("data/pgs/", trait, "_regenie_1_scores.sscore"))
  temp$pgs <- temp$SCORE1_AVG
  temp <- temp %>% select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG))
  temp$phenotype <- rep(trait, nrow(temp))
  pgs_df_regenie <- rbind.data.frame(pgs_df_regenie, temp)
  rm(temp)
}

pgs_df_regenie <- left_join(group_non_pgs_df, pgs_df_regenie, by = c("#FID", "IID", "phenotype"))
pgs_df_regenie <- pgs_df_regenie %>% select(`#FID`, IID, sex, age, age_sq, age_sex, age_sq_sex, array_type,
                                                      pc_dist, pgs, phenotype, phenotype_value, weighted_pc_groups, group_close_to_gwas)

pgs_df_regenie <- pgs_df_regenie %>% mutate(
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
temp <- pgs_df_regenie %>% filter(phenotype == pheno[1] & weighted_pc_groups == 1)
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
for(i in 2:250){
  temp2 <- pgs_df_regenie %>% filter(phenotype == pheno[1] & weighted_pc_groups == i)
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
  for(i in 1:250){
    temp2 <- pgs_df_regenie %>% filter(phenotype == trait & weighted_pc_groups == i)
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

pgs_df_regenie <- temp %>% arrange(IID)

# X = Y - Y_hat
pgs_df_regenie$X <- pgs_df_regenie$phenotype_value - pgs_df_regenie$Y_hat

# Do X ~ gen_dist_polynomial + sex + sex_gen_dist across all bins
X_mod <- pgs_df_regenie %>%
  group_by(phenotype) %>%
  do(lm_X = lm(X ~ pc_dist + pc_dist_sq + pc_dist_3 + pc_dist_4 + pc_dist_5 + pc_dist_6 +
                 pc_dist_7 + pc_dist_8 + pc_dist_9 + pc_dist_10 + pc_dist_11 + pc_dist_12 +
                 pc_dist_13 + pc_dist_14 + pc_dist_15 + pc_dist_16 + pc_dist_17 + pc_dist_18 +
                 pc_dist_19 + pc_dist_20 + sex + sex_pc_dist, data = .))

# Apply this model
temp <- pgs_df_regenie %>% filter(phenotype == pheno[1])
temp$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_regenie %>% filter(phenotype == trait)
  temp2$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_regenie <- temp %>% arrange(IID)

# Z = X - X_hat
pgs_df_regenie$Z <- pgs_df_regenie$X - pgs_df_regenie$X_hat

# Z ~ PGS
Z_mod <- pgs_df_regenie %>%
  group_by(phenotype) %>%
  do(lm_Z = lm(Z ~ pgs, data = .))

# Apply this model
temp <- pgs_df_regenie %>% filter(phenotype == pheno[1])
temp$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_regenie %>% filter(phenotype == trait)
  temp2$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_regenie <- temp %>% arrange(IID)

# Prediction error = (Z - Z_hat) ^ 2
pgs_df_regenie$pred_error = (pgs_df_regenie$Z - pgs_df_regenie$Z_hat) ^ 2

pgs_df_regenie <- pgs_df_regenie %>% arrange(IID, phenotype)

pgs_df_regenie %>% write_tsv("data/pgs_pred/ind_pgs_df_regenie.tsv")


# PGS with PRS-CS
# Read in PGS for threshold 1 (p < 1e-5)
pgs_df_prscs <- data.frame()
for(trait in pheno){
  temp <- read_tsv(paste0("data/pgs/", trait, "_prscs_scores.sscore"))
  temp$pgs <- temp$SCORE1_AVG
  temp <- temp %>% select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG))
  temp$phenotype <- rep(trait, nrow(temp))
  pgs_df_prscs <- rbind.data.frame(pgs_df_prscs, temp)
  rm(temp)
}

pgs_df_prscs <- left_join(group_non_pgs_df, pgs_df_prscs, by = c("#FID", "IID", "phenotype"))
pgs_df_prscs <- pgs_df_prscs %>% select(`#FID`, IID, sex, age, age_sq, age_sex, age_sq_sex, array_type,
                                            pc_dist, pgs, phenotype, phenotype_value, weighted_pc_groups, group_close_to_gwas)

pgs_df_prscs <- pgs_df_prscs %>% mutate(
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
temp <- pgs_df_prscs %>% filter(phenotype == pheno[1] & weighted_pc_groups == 1)
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
for(i in 2:250){
  temp2 <- pgs_df_prscs %>% filter(phenotype == pheno[1] & weighted_pc_groups == i)
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
  for(i in 1:250){
    temp2 <- pgs_df_prscs %>% filter(phenotype == trait & weighted_pc_groups == i)
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

pgs_df_prscs <- temp %>% arrange(IID)

# X = Y - Y_hat
pgs_df_prscs$X <- pgs_df_prscs$phenotype_value - pgs_df_prscs$Y_hat

# Do X ~ gen_dist_polynomial + sex + sex_gen_dist across all bins
X_mod <- pgs_df_prscs %>%
  group_by(phenotype) %>%
  do(lm_X = lm(X ~ pc_dist + pc_dist_sq + pc_dist_3 + pc_dist_4 + pc_dist_5 + pc_dist_6 +
                 pc_dist_7 + pc_dist_8 + pc_dist_9 + pc_dist_10 + pc_dist_11 + pc_dist_12 +
                 pc_dist_13 + pc_dist_14 + pc_dist_15 + pc_dist_16 + pc_dist_17 + pc_dist_18 +
                 pc_dist_19 + pc_dist_20 + sex + sex_pc_dist, data = .))

# Apply this model
temp <- pgs_df_prscs %>% filter(phenotype == pheno[1])
temp$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_prscs %>% filter(phenotype == trait)
  temp2$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_prscs <- temp %>% arrange(IID)

# Z = X - X_hat
pgs_df_prscs$Z <- pgs_df_prscs$X - pgs_df_prscs$X_hat

# Z ~ PGS
Z_mod <- pgs_df_prscs %>%
  group_by(phenotype) %>%
  do(lm_Z = lm(Z ~ pgs, data = .))

# Apply this model
temp <- pgs_df_prscs %>% filter(phenotype == pheno[1])
temp$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_prscs %>% filter(phenotype == trait)
  temp2$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_prscs <- temp %>% arrange(IID)

# Prediction error = (Z - Z_hat) ^ 2
pgs_df_prscs$pred_error = (pgs_df_prscs$Z - pgs_df_prscs$Z_hat) ^ 2

pgs_df_prscs <- pgs_df_prscs %>% arrange(IID, phenotype)

pgs_df_prscs %>% write_tsv("data/pgs_pred/ind_pgs_df_prscs.tsv")


# PC distance with 16 PCs
# Read in group level file for covariates (16 PCs)
group_non_pgs_df_16 <- read_tsv("data/pgs_pred/group_non_pgs_df_16.tsv")

# Read in PGS for threshold 1 (p < 1e-5)
pgs_df_16 <- data.frame()
for(trait in pheno){
  temp <- read_tsv(paste0("data/pgs/", trait, "_1_scores.sscore"))
  temp$pgs <- temp$SCORE1_AVG
  temp <- temp %>% select(-c(ALLELE_CT, NAMED_ALLELE_DOSAGE_SUM, SCORE1_AVG))
  temp$phenotype <- rep(trait, nrow(temp))
  pgs_df_16 <- rbind.data.frame(pgs_df_16, temp)
  rm(temp)
}

pgs_df_16 <- left_join(group_non_pgs_df_16, pgs_df_16, by = c("#FID", "IID", "phenotype"))
pgs_df_16 <- pgs_df_16 %>% select(`#FID`, IID, sex, age, age_sq, age_sex, age_sq_sex, array_type,
                            pc_dist, pgs, phenotype, phenotype_value, weighted_pc_groups, group_close_to_gwas)

pgs_df_16 <- pgs_df_16 %>% mutate(
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
temp <- pgs_df_16 %>% filter(phenotype == pheno[1] & weighted_pc_groups == 1)
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
for(i in 2:250){
  temp2 <- pgs_df_16 %>% filter(phenotype == pheno[1] & weighted_pc_groups == i)
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
  for(i in 1:250){
    temp2 <- pgs_df_16 %>% filter(phenotype == trait & weighted_pc_groups == i)
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

pgs_df_16 <- temp %>% arrange(IID)

# X = Y - Y_hat
pgs_df_16$X <- pgs_df_16$phenotype_value - pgs_df_16$Y_hat

# Do X ~ gen_dist_polynomial + sex + sex_gen_dist across all bins
X_mod <- pgs_df_16 %>%
  group_by(phenotype) %>%
  do(lm_X = lm(X ~ pc_dist + pc_dist_sq + pc_dist_3 + pc_dist_4 + pc_dist_5 + pc_dist_6 +
                 pc_dist_7 + pc_dist_8 + pc_dist_9 + pc_dist_10 + pc_dist_11 + pc_dist_12 +
                 pc_dist_13 + pc_dist_14 + pc_dist_15 + pc_dist_16 + pc_dist_17 + pc_dist_18 +
                 pc_dist_19 + pc_dist_20 + sex + sex_pc_dist, data = .))

# Apply this model
temp <- pgs_df_16 %>% filter(phenotype == pheno[1])
temp$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_16 %>% filter(phenotype == trait)
  temp2$X_hat <- unname(predict(X_mod$lm_X[X_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_16 <- temp %>% arrange(IID)

# Z = X - X_hat
pgs_df_16$Z <- pgs_df_16$X - pgs_df_16$X_hat

# Z ~ PGS
Z_mod <- pgs_df_16 %>%
  group_by(phenotype) %>%
  do(lm_Z = lm(Z ~ pgs, data = .))

# Apply this model
temp <- pgs_df_16 %>% filter(phenotype == pheno[1])
temp$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == pheno[1]][[1]], newdata = temp))
for(trait in pheno[2:15]){
  temp2 <- pgs_df_16 %>% filter(phenotype == trait)
  temp2$Z_hat <- unname(predict(Z_mod$lm_Z[Z_mod$phenotype == trait][[1]], newdata = temp2))
  temp <- rbind.data.frame(temp, temp2)
}
pgs_df_16 <- temp %>% arrange(IID)

# Prediction error = (Z - Z_hat) ^ 2
pgs_df_16$pred_error = (pgs_df_16$Z - pgs_df_16$Z_hat) ^ 2

pgs_df_16 <- pgs_df_16 %>% arrange(IID, phenotype)

pgs_df_16 %>% write_tsv("data/pgs_pred/ind_pgs_df_16.tsv")
