library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(stringr)

# Use the file for the for the ones that don't require non_pgs_df to be changed
non_pgs_df <- read_tsv("data/pgs_pred/group_non_pgs_df.tsv")

partial.R2.helper <- function(nested.lm, ref.lm){
  a <- anova(nested.lm)
  b <- anova(ref.lm)
  length.ref <- length(attributes(ref.lm$terms)$"dataClasses")
  length.nested <- length(attributes(nested.lm$terms)$"dataClasses")
  if(length.nested > length.ref) stop("Specify nested model first in arguements")
  if(length.ref - length.nested > 1) stop("Reference and nested model should only differ with repsect to the presence/absence of one predictor")
  SSE.wo <- tail(a$"Sum Sq", 1)
  SSE.with <- tail(b$"Sum Sq", 1)
  P.R2<-(SSE.wo-SSE.with)/SSE.wo
  return(P.R2)
}

partial.R2 <- function(nested.lm.list, ref.lm.list){
  res = c()
  for(i in 1:length(nested.lm.list)){
    res = c(res, partial.R2.helper(nested.lm.list[[i]], ref.lm.list[[i]]))
  }
  return(res)
}

load_non_pgs_df <- function(num_bins) {
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar.tsv')

  # Drop all PC columns
  covar_df <- covar_df %>% select(-starts_with("PC"))
  
  # Load array type
  array_df <- read_table("data/extracted_data_fields/array_type.txt")
  colnames(array_df) = c("IID", "array_batch")
  array_df$array_type = ifelse(array_df$array_batch < 0, "BiLEVE",
                               ifelse(array_df$array_batch > 0 & array_df$array_batch < 100, "Axiom",
                                      ifelse(array_df$array_batch == 1000, "BiLEVE",
                                             ifelse(array_df$array_batch == 2000, "Axiom", array_df$array_batch))))
  
  covar_df <- covar_df %>% filter(pop != "WB_GWAS")
  
  phenotypes_df <- read_tsv('data/phenotypes/phenotypes.tsv')
  
  # Load the individuals 
  population_files <- c('data/ukb_populations/nwb_all_id.txt', 
                        'data/ukb_populations/wb_gwas_id.txt')
  populations_df <- data.frame()
  for (file in population_files) {
    pop_df <- read_delim(file, delim = ' ', trim_ws = T,
                         col_types = c('#FID' = col_integer(), 'IID' = col_integer())) %>%
      mutate(population = str_extract(file, '(?<=data/ukb_populations/)[a-z]+_[a-z]+'))
    populations_df <- bind_rows(populations_df, pop_df)
  }
  pc_values <- read_tsv('data/pca/pc_dist_best_pred_std.tsv')
  pc_values <- pc_values[, c(2:3)]
  
  # Combine all the above tables into a table that will be joined with PGS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(array_df, by = "IID") %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    left_join(pc_values, by = c("IID")) %>%
    pivot_longer(Height:LDL, names_to = 'phenotype', values_to = 'phenotype_value')
  
  # Divides into equally-sized components
  phenotypes_df <- phenotypes_df %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  return(phenotypes_df)
}

# For stratifying individuals by Townsend index. num_strata can only be either 2 or 5
load_non_pgs_df_townsend <- function(num_bins, num_strata) {
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar.tsv')

  # Drop all PC columns
  covar_df <- covar_df %>% select(-starts_with("PC"))
  
  # Load array type
  array_df <- read_table("data/extracted_data_fields/array_type.txt")
  colnames(array_df) = c("IID", "array_batch")
  array_df$array_type = ifelse(array_df$array_batch < 0, "BiLEVE",
                               ifelse(array_df$array_batch > 0 & array_df$array_batch < 100, "Axiom",
                                      ifelse(array_df$array_batch == 1000, "BiLEVE",
                                             ifelse(array_df$array_batch == 2000, "Axiom", array_df$array_batch))))
  
  covar_df <- covar_df %>% filter(pop != "WB_GWAS")
  
  phenotypes_df <- read_tsv('data/phenotypes/phenotypes.tsv')
  
  # Load the individuals 
  population_files <- c('data/ukb_populations/nwb_all_id.txt', 
                        'data/ukb_populations/wb_gwas_id.txt')
  populations_df <- data.frame()
  for (file in population_files) {
    pop_df <- read_delim(file, delim = ' ', trim_ws = T,
                         col_types = c('#FID' = col_integer(), 'IID' = col_integer())) %>%
      mutate(population = str_extract(file, '(?<=data/ukb_populations/)[a-z]+_[a-z]+'))
    populations_df <- bind_rows(populations_df, pop_df)
  }
  pc_values <- read_tsv('data/pca/pc_dist_best_pred_std.tsv')
  pc_values <- pc_values[, c(2:3)]

  townsend <- read_table("data/extracted_phenotypes/townsend.txt")
  colnames(townsend) <- c("IID", "townsend")
  townsend <- townsend %>% filter(IID %in% covar_df$IID)

  townsend$townsend_rank <- ntile(townsend$townsend, num_strata)
  
  # Combine all the above tables into a table that will be joined with PGS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(array_df, by = "IID") %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    left_join(townsend, by = c('IID')) %>%
    left_join(pc_values, by = c("IID")) %>%
    pivot_longer(Height:LDL, names_to = 'phenotype', values_to = 'phenotype_value')
  
  if(num_strata == 2){
    phenotypes_df_1 <- phenotypes_df %>% filter(townsend_rank == 1)
    phenotypes_df_2 <- phenotypes_df %>% filter(townsend_rank == 2)
  
  phenotypes_df_1 <- phenotypes_df_1 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df_2 <- phenotypes_df_2 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df <- rbind.data.frame(phenotypes_df_1, phenotypes_df_2)
    
  } else if(num_strata == 5){
    phenotypes_df_1 <- phenotypes_df %>% filter(townsend_rank == 1)
    phenotypes_df_2 <- phenotypes_df %>% filter(townsend_rank == 2)
    phenotypes_df_3 <- phenotypes_df %>% filter(townsend_rank == 3)
    phenotypes_df_4 <- phenotypes_df %>% filter(townsend_rank == 4)
    phenotypes_df_5 <- phenotypes_df %>% filter(townsend_rank == 5)
  
  phenotypes_df_1 <- phenotypes_df_1 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df_2 <- phenotypes_df_2 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df_3 <- phenotypes_df_3 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df_4 <- phenotypes_df_4 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df_5 <- phenotypes_df_5 %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  phenotypes_df <- rbind.data.frame(phenotypes_df_1, phenotypes_df_2, phenotypes_df_3, 
                                    phenotypes_df_4, phenotypes_df_5)
  
  }

  pc_groups <- phenotypes_df %>% select(weighted_pc_groups, median_pc_dist) %>% unique() %>%
    arrange(median_pc_dist)
  pc_groups$weighted_pc_groups_new <- 1:(num_bins * num_strata)
  
  phenotypes_df <- phenotypes_df %>%
    left_join(pc_groups, by = c("weighted_pc_groups", "median_pc_dist")) %>%
    select(-weighted_pc_groups) %>%
    rename(weighted_pc_groups = weighted_pc_groups_new)
  
  return(phenotypes_df)
}

get_r2_values <- function(df, group_var) {
  df$group <- df %>% select(any_of(group_var)) %>% as.data.frame()
  
  df1 = data.frame()
  for(g in unique(df$group$weighted_pc_groups)){
    for(p in unique(df$phenotype)){
      for(t in unique(df$threshold)){
        temp <- df %>% filter(group == g & phenotype == p & threshold == t)
        meta <- as.data.frame(t(c(g, p, t)))
        colnames(meta) <- c("weighted_pc_groups", "phenotype", "threshold")
        meta$weighted_pc_groups <- as.numeric(meta$weighted_pc_groups)
        meta$threshold <- as.numeric(meta$threshold)
        
        # Use the array type as a covariate only when both types have more than 2 samples
        if(length(temp$array_type[temp$array_type == "Axiom"]) > 2 & length(temp$array_type[temp$array_type == "BiLEVE"]) > 2){
          meta_2 <- meta %>%
            do(nested = lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex + array_type, data = temp),
               full = lm(phenotype_value ~ pgs + sex + age + age_sq + age_sex + age_sq_sex + array_type, data = temp))
        } else{
          meta_2 <- meta %>%
            do(nested = lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex, data = temp),
               full = lm(phenotype_value ~ pgs + sex + age + age_sq + age_sex + age_sq_sex, data = temp))
        }
        df1 <- rbind.data.frame(df1, cbind.data.frame(meta, meta_2))
      }
    }
  }
  
  df1 <- df1 %>%
    mutate(partial = partial.R2(nested, full)) %>%
    select(-nested, -full)
  df1$group <- 1:nrow(df1)
  colnames(df1)[1] <- group_var
  return(df1)
}

# Calculate PGS prediction accuracy
make_pgs_evaluation_df <- function(non_pgs_df, group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  pgs_df <- data.frame()
  
  files <- list.files(path = 'data/pgs', pattern = '[a-zA-Z]+_[0-9]_scores.sscore', 
                      full.names = TRUE)
  files <- files[!grepl("disease|array|regenie|prscs", files)]
  
  for (file in files) {
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', pgs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z_]+(?=_)'),
        threshold = str_extract_all(string = file, pattern = '[0-4]')[[1]] %>% as.integer
      ) %>%
      inner_join(non_pgs_df, by =  c('#FID', 'IID','phenotype')) %>%
      get_r2_values(., group_var = group_var_as_string)
    
    pgs_df <- bind_rows(pgs_df, this_df)
  }
  return(pgs_df)
}

# Data frame without PGS information
non_pgs_df_500_bins <- load_non_pgs_df(num_bins = 500)
non_pgs_df_500_bins <- non_pgs_df_500_bins %>% filter(!is.na(weighted_pc_groups))
non_pgs_df_500_bins <- non_pgs_df_500_bins %>%
  mutate(array_type = as.factor(array_type))

# Data frame with PGS information
pgs_df_500_bins  <- make_pgs_evaluation_df(non_pgs_df_500_bins, "weighted_pc_groups")
pgs_df_500_bins <- pgs_df_500_bins %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)

# Data frame without PGS information
# 1/2 of 250 bins because of 2 stata
non_pgs_df_townsend_2_strata <- load_non_pgs_df_townsend(num_bins = 125, num_strata = 2)
non_pgs_df_townsend_2_strata <- non_pgs_df_townsend_2_strata %>% filter(!is.na(weighted_pc_groups))
non_pgs_df_townsend_2_strata <- non_pgs_df_townsend_2_strata %>%
  mutate(array_type = as.factor(array_type))

# Data frame with PGS information
pgs_df_townsend_2_strata  <- make_pgs_evaluation_df(non_pgs_df_townsend_2_strata, "weighted_pc_groups")
pgs_df_townsend_2_strata <- pgs_df_townsend_2_strata %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)

# 1/5 of 250 bins because of 5 stata
non_pgs_df_townsend_5_strata <- load_non_pgs_df_townsend(num_bins = 50, num_strata = 5)
non_pgs_df_townsend_5_strata <- non_pgs_df_townsend_5_strata %>% filter(!is.na(weighted_pc_groups))
non_pgs_df_townsend_5_strata <- non_pgs_df_townsend_5_strata %>%
  mutate(array_type = as.factor(array_type))

# Data frame with PGS information
pgs_df_townsend_5_strata  <- make_pgs_evaluation_df(non_pgs_df_townsend_5_strata, "weighted_pc_groups")
pgs_df_townsend_5_strata <- pgs_df_townsend_5_strata %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)


# GWAS with assessment center and genotype array
# Calculate PGS prediction accuracy
make_pgs_evaluation_df_array_center <- function(non_pgs_df, group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  pgs_df <- data.frame()
  
  for (file in list.files(path = 'data/pgs',
                         pattern = '[a-zA-Z]+_array_center_[0-9]_scores.sscore', full.names = T)){
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', pgs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z_]+(?=_array_center)'),
        threshold = str_extract_all(string = file, pattern = '[0-4]')[[1]] %>% as.integer
      ) %>%
      inner_join(non_pgs_df, by =  c('#FID', 'IID','phenotype')) %>%
      get_r2_values(., group_var = group_var_as_string)
    
    pgs_df <- bind_rows(pgs_df, this_df)
  }
  return(pgs_df)
}

# Data frame with PGS information
pgs_df_array_center  <- make_pgs_evaluation_df_array_center(non_pgs_df, "weighted_pc_groups")
pgs_df_array_center <- pgs_df_array_center %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)


# GWAS with regenie
# Calculate PGS prediction accuracy
make_pgs_evaluation_df_regenie <- function(non_pgs_df, group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  pgs_df <- data.frame()
  
  for (file in list.files(path = 'data/pgs',
                          pattern = '[a-zA-Z]+_regenie_[0-9]_scores.sscore', full.names = T)){
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', pgs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z_]+(?=_regenie)'),
        threshold = str_extract_all(string = file, pattern = '[0-4]')[[1]] %>% as.integer
      ) %>%
      inner_join(non_pgs_df, by =  c('#FID', 'IID','phenotype')) %>%
      get_r2_values(., group_var = group_var_as_string)
    
    pgs_df <- bind_rows(pgs_df, this_df)
  }
  return(pgs_df)
}

# Data frame with PGS information
pgs_df_regenie  <- make_pgs_evaluation_df_regenie(non_pgs_df, "weighted_pc_groups")
pgs_df_regenie <- pgs_df_regenie %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)

# PGS with PRS-CS
get_r2_values_prscs <- function(df, group_var) {
  df$group <- df %>% select(any_of(group_var)) %>% as.data.frame()
  
  df1 = data.frame()
  for(g in unique(df$group$weighted_pc_groups)){
    for(p in unique(df$phenotype)){
      temp <- df %>% filter(group == g & phenotype == p)
      meta <- as.data.frame(t(c(g, p)))
      colnames(meta) <- c("weighted_pc_groups", "phenotype")
      meta$weighted_pc_groups <- as.numeric(meta$weighted_pc_groups)
      
      # Use the array type as a covariate only when both types have more than 2 samples
      if(length(temp$array_type[temp$array_type == "Axiom"]) > 2 & length(temp$array_type[temp$array_type == "BiLEVE"]) > 2){
        meta_2 <- meta %>%
          do(nested = lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex + array_type, data = temp),
             full = lm(phenotype_value ~ pgs + sex + age + age_sq + age_sex + age_sq_sex + array_type, data = temp))
      } else{
        meta_2 <- meta %>%
          do(nested = lm(phenotype_value ~ sex + age + age_sq + age_sex + age_sq_sex, data = temp),
             full = lm(phenotype_value ~ pgs + sex + age + age_sq + age_sex + age_sq_sex, data = temp))
      }
      df1 <- rbind.data.frame(df1, cbind.data.frame(meta, meta_2))
    }
  }
  
  df1 <- df1 %>%
    mutate(partial = partial.R2(nested, full)) %>%
    select(-nested, -full)
  df1$group <- 1:nrow(df1)
  colnames(df1)[1] <- group_var
  return(df1)
}

# Calculate PGS prediction accuracy
make_pgs_evaluation_df_prscs <- function(non_pgs_df, group_var_as_string) {
  # Combine partial R^2 values for each phenotype and threshold
  col_types = c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  pgs_df <- data.frame()
  
  for (file in list.files(path = 'data/pgs',
                          pattern = '[a-zA-Z]+_prscs_scores.sscore', full.names = T)){
    print(file)
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', pgs = 'SCORE1_AVG') %>%
      mutate(
        phenotype = str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z_]+(?=_prscs)')
      ) %>%
      inner_join(non_pgs_df, by =  c('#FID', 'IID','phenotype')) %>%
      get_r2_values_prscs(., group_var = group_var_as_string)
    
    pgs_df <- bind_rows(pgs_df, this_df)
  }
  return(pgs_df)
}

# Data frame with PGS information
pgs_df_prscs <- make_pgs_evaluation_df_prscs(non_pgs_df, "weighted_pc_groups")
pgs_df_prscs <- pgs_df_prscs %>% 
  arrange(phenotype, weighted_pc_groups) %>%
  mutate(group_number = weighted_pc_groups)


# Genetic distance with 16 PCs
load_non_pgs_df_16 <- function(num_bins) {
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar.tsv')

  # Drop all PC columns
  covar_df <- covar_df %>% select(-starts_with("PC"))
  
  # Load array type
  array_df <- read_table("data/extracted_data_fields/array_type.txt")
  colnames(array_df) = c("IID", "array_batch")
  array_df$array_type = ifelse(array_df$array_batch < 0, "BiLEVE",
                               ifelse(array_df$array_batch > 0 & array_df$array_batch < 100, "Axiom",
                                      ifelse(array_df$array_batch == 1000, "BiLEVE",
                                             ifelse(array_df$array_batch == 2000, "Axiom", array_df$array_batch))))
  
  covar_df <- covar_df %>% filter(pop != "WB_GWAS")
  
  phenotypes_df <- read_tsv('data/phenotypes/phenotypes.tsv')
  
  # Load the individuals 
  population_files <- c('data/ukb_populations/nwb_all_id.txt', 
                        'data/ukb_populations/wb_gwas_id.txt')
  populations_df <- data.frame()
  for (file in population_files) {
    pop_df <- read_delim(file, delim = ' ', trim_ws = T,
                         col_types = c('#FID' = col_integer(), 'IID' = col_integer())) %>%
      mutate(population = str_extract(file, '(?<=data/ukb_populations/)[a-z]+_[a-z]+'))
    populations_df <- bind_rows(populations_df, pop_df)
  }
  pc_values <- read_tsv('data/pca/pc_dist_16.tsv')
  pc_values <- pc_values[, c(2:3)]
  
  # Combine all the above tables into a table that will be joined with PGS information for
  # phenotype-threshold combinations.
  phenotypes_df <- phenotypes_df %>%
    inner_join(covar_df, by = c('#FID', 'IID')) %>%
    inner_join(array_df, by = "IID") %>%
    inner_join(populations_df, by = c('#FID', 'IID')) %>%
    left_join(pc_values, by = c("IID")) %>%
    pivot_longer(Height:LDL, names_to = 'phenotype', values_to = 'phenotype_value')
  
  # Divides into equally-sized components
  phenotypes_df <- phenotypes_df %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  return(phenotypes_df)
}

# Data frame without PGS information
non_pgs_df_16 <- load_non_pgs_df_16(num_bins = 250)
non_pgs_df_16 <- non_pgs_df_16 %>% filter(!is.na(weighted_pc_groups))
non_pgs_df_16 <- non_pgs_df_16 %>%
  mutate(array_type = as.factor(array_type))

# Data frame with PGS information
pgs_df_16  <- make_pgs_evaluation_df(non_pgs_df_16, "weighted_pc_groups")
pgs_df_16 <- pgs_df_16 %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)



# Calculate the median PC distance for each group
get_median_pc <- function(file){
  median_pc_values <- file %>% 
    filter(phenotype =="Height") %>%
    select(weighted_pc_groups, pc_dist) %>%
    group_by(weighted_pc_groups) %>%
    dplyr::summarize(median_pc = median(pc_dist),
                     count = n()) %>%
    na.omit() %>%
    as.data.frame()
  return(median_pc_values)
}

# Order the bins by how close to genetic distance = 1 (mean distance to the GWAS centroid of the GWAS group)
median_pc_500_bins <- get_median_pc(non_pgs_df_500_bins)
median_pc_500_bins$group_close_to_gwas <- abs(median_pc_500_bins$median_pc - 1)
median_pc_500_bins <- median_pc_500_bins %>% arrange(group_close_to_gwas)
median_pc_500_bins$group_close_to_gwas <- 1:500
median_pc_500_bins <- median_pc_500_bins %>% arrange(weighted_pc_groups)
non_pgs_df_500_bins <- non_pgs_df_500_bins %>% left_join(median_pc_500_bins[, c(1:2, 4)], by = "weighted_pc_groups")
pgs_df_500_bins <- pgs_df_500_bins %>% left_join(median_pc_500_bins[, c(1:2, 4)], by = "weighted_pc_groups")

median_pc <- get_median_pc(non_pgs_df)
median_pc$group_close_to_gwas <- abs(median_pc$median_pc - 1)
median_pc <- median_pc %>% arrange(group_close_to_gwas)
median_pc$group_close_to_gwas <- 1:250
median_pc <- median_pc %>% arrange(weighted_pc_groups)
non_pgs_df <- non_pgs_df %>% left_join(median_pc[, c(1:2, 4)], by = "weighted_pc_groups")

pgs_df_array_center <- pgs_df_array_center %>% left_join(median_pc[, c(1:2, 4)], by = "weighted_pc_groups")

pgs_df_regenie <- pgs_df_regenie %>% left_join(median_pc[, c(1:2, 4)], by = "weighted_pc_groups")

pgs_df_prscs <- pgs_df_prscs %>% left_join(median_pc[, c(1:2, 4)], by = "weighted_pc_groups")

median_pc_16 <- get_median_pc(non_pgs_df_16)
median_pc_16$group_close_to_gwas <- abs(median_pc_16$median_pc - 1)
median_pc_16 <- median_pc_16 %>% arrange(group_close_to_gwas)
median_pc_16$group_close_to_gwas <- 1:250
median_pc_16 <- median_pc_16 %>% arrange(weighted_pc_groups)
non_pgs_df_16 <- non_pgs_df_16 %>% left_join(median_pc_16[, c(1:2, 4)], by = "weighted_pc_groups")
pgs_df_16 <- pgs_df_16 %>% left_join(median_pc_16[, c(1:2, 4)], by = "weighted_pc_groups")

# Save the files
non_pgs_df_500_bins %>% write_tsv("data/pgs_pred/group_non_pgs_df_500_bins.tsv")
pgs_df_500_bins %>% write_tsv("data/pgs_pred/group_pgs_df_500_bins.tsv")

non_pgs_df_townsend_2_strata %>% write_tsv("data/pgs_pred/group_non_pgs_df_townsend_2_strata.tsv")
pgs_df_townsend_2_strata %>% write_tsv("data/pgs_pred/group_pgs_df_townsend_2_strata.tsv")

non_pgs_df_townsend_5_strata %>% write_tsv("data/pgs_pred/group_non_pgs_df_townsend_5_strata.tsv")
pgs_df_townsend_5_strata %>% write_tsv("data/pgs_pred/group_pgs_df_townsend_5_strata.tsv")

pgs_df_array_center %>% write_tsv("data/pgs_pred/group_pgs_df_array_center.tsv")

pgs_df_regenie %>% write_tsv("data/pgs_pred/group_pgs_df_regenie.tsv")

pgs_df_prscs %>% write_tsv("data/pgs_pred/group_pgs_df_prscs.tsv")

non_pgs_df_16 %>% write_tsv("data/pgs_pred/group_non_pgs_df_16.tsv")
pgs_df_16 %>% write_tsv("data/pgs_pred/group_pgs_df_16.tsv")
