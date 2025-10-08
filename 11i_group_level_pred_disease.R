library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(stringr)

get_precision_recall <- function(df, phenotype, thres){
  
  df <- df %>% mutate(phenotype_value_pred = ifelse(prs >= thres, 1, 0))
  
  # For Alzheimer's, make 25 bins instead of 1,000 due to the low number of
  # positive individuals
  if(phenotype == "Alzheimer"){
    df1 <- df %>% mutate(weighted_pc_groups = pc_dist %>% ntile(25)) %>%
      filter(!is.na(phenotype_value), !is.na(phenotype_value_pred)) %>%
      group_by(phenotype, threshold, weighted_pc_groups) %>%
      summarize(TP = sum(phenotype_value == 1 & phenotype_value_pred == 1),
                FP = sum(phenotype_value_pred == 1 & phenotype_value == 0),
                FN = sum(phenotype_value_pred == 0 & phenotype_value == 1),
                precision = TP / (TP + FP),
                recall = TP / (TP + FN),
                f1_score = 2 * (precision * recall) / (precision + recall),
                median_pc_dist = median(pc_dist))
  } else{
    df1 <- df %>% mutate(new_pc_groups = weighted_pc_groups) %>%
      filter(!is.na(phenotype_value), !is.na(phenotype_value_pred)) %>%
      group_by(phenotype, threshold, weighted_pc_groups) %>%
      summarize(TP = sum(phenotype_value == 1 & phenotype_value_pred == 1),
                FP = sum(phenotype_value_pred == 1 & phenotype_value == 0),
                FN = sum(phenotype_value_pred == 0 & phenotype_value == 1),
                precision = TP / (TP + FP),
                recall = TP / (TP + FN),
                f1_score = 2 * (precision * recall) / (precision + recall),
                median_pc_dist = median(pc_dist))
  }
  
  
  df1 <- df1 %>%
    select(-c(TP, FP, FN, f1_score, median_pc_dist))
  
  return(df1)
}

load_non_pgs_df_disease <- function(num_bins) {
  # Load covariates (age, sex, age-sex interactions, PC1, ..., PC20)
  covar_df <- read_tsv('data/ukb_merged/covar.tsv')

  # Drop the PC columns
  covar_df <- covar_df %>% select(-starts_with("PC"))
  
  # Load array type
  array_df <- read_table("data/extracted_data_fields/array_type.txt")
  colnames(array_df) = c("IID", "array_batch")
  array_df$array_type = ifelse(array_df$array_batch < 0, "BiLEVE",
                               ifelse(array_df$array_batch > 0 & array_df$array_batch < 100, "Axiom",
                                      ifelse(array_df$array_batch == 1000, "BiLEVE",
                                             ifelse(array_df$array_batch == 2000, "Axiom", array_df$array_batch))))
  
  covar_df <- covar_df %>% filter(pop != "WB_GWAS")
  
  phenotypes_df <- read_tsv('data/phenotypes/phenotypes_disease.tsv')
  
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
    pivot_longer(Alzheimer:Asthma, names_to = 'phenotype', values_to = 'phenotype_value')
  
  # Divides into equally-sized components
  phenotypes_df <- phenotypes_df %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(num_bins))
  
  return(phenotypes_df)
}

# Calculate precision and recall
get_precision_recall_disease <- function(non_prs_df) {
  col_types <- c('#FID' = col_integer(), 'IID' = col_integer(), 'ALLELE_CT' = col_integer(),
                 'NAMED_ALLELE_DOSAGE_SUM' = col_double(), 'SCORE1_AVG' = col_double())
  prs_df <- data.frame()
  
  for (file in list.files(path = 'data/pgs',
                          pattern = '[a-zA-Z]+_disease_[0-9]_scores.sscore', full.names = T)) {
    print(file)
    
    this_df <- read_tsv(file, col_types = col_types) %>%
      filter(IID > 0) %>%
      select('#FID', 'IID', prs = 'SCORE1_AVG') 
    
    gwas_prs <- read_tsv(gsub("_disease_", "_disease_gwas_", file))
    
    # Got these thresholds by maximizing the F1 scores in the GWAS groups for each trait
    perc = ifelse(str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z0-9]+(?=_)') == "Alzheimer", 0.95,
                  ifelse(str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z0-9]+(?=_)') == "T2D", 0.75, 
                         ifelse(str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z0-9]+(?=_)') == "Asthma", 0.65, NA)))
    thres <- quantile(gwas_prs$SCORE1_AVG, perc, na.rm = TRUE)
    
    this_df <- this_df %>% mutate(
      phenotype = str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z0-9]+(?=_)'),
      threshold = ifelse(phenotype == "T2D", str_extract_all(string = file, pattern = '[0-4]')[[1]][2] %>% as.integer, 
                         str_extract_all(string = file, pattern = '[0-4]')[[1]] %>% as.integer)
    ) %>%
      inner_join(non_prs_df, by = c('#FID', 'IID', 'phenotype')) %>%
      get_precision_recall(., str_extract(string = file, pattern = '(?<=data/pgs/)[A-Za-z0-9]+(?=_)'), thres)
    
    prs_df <- bind_rows(prs_df, this_df)
  }
  
  return(prs_df)
}

# Data frame without PGS information
non_pgs_df <- load_non_pgs_df_disease(num_bins = 250)
non_pgs_df <- non_pgs_df %>%
  mutate(array_type = as.factor(array_type))

# Data frame with PGS information
pgs_df  <- get_precision_recall_disease(non_pgs_df)
pgs_df <- pgs_df %>% 
  arrange(phenotype, weighted_pc_groups, threshold) %>%
  mutate(group_number = weighted_pc_groups)

# Calculate the median PC distance for each group
get_median_pc_disease = function(file){
  median_pc_values <- file %>% 
    filter(phenotype =="Alzheimer") %>%
    mutate(weighted_pc_groups = pc_dist %>% ntile(25)) %>%
    select(weighted_pc_groups, pc_dist) %>%
    group_by(weighted_pc_groups) %>%
    dplyr::summarize(median_pc = median(pc_dist),
                     count = n()) %>%
    na.omit() %>%
    as.data.frame()
  return(median_pc_values)
}

get_median_pc_disease_2 = function(file){
  median_pc_values <- file %>% 
    filter(phenotype =="T2D") %>%
    select(weighted_pc_groups, pc_dist) %>%
    group_by(weighted_pc_groups) %>%
    dplyr::summarize(median_pc = median(pc_dist),
                     count = n()) %>%
    na.omit() %>%
    as.data.frame()
  return(median_pc_values)
}

# Order the bins by how close to genetic distance = 1 (mean distance to the GWAS centroid of the GWAS group)
median_pc_alzheimer <- get_median_pc_disease(non_pgs_df)
median_pc_alzheimer$group_close_to_gwas <- abs(median_pc_alzheimer$median_pc - 1)
median_pc_alzheimer <- median_pc_alzheimer %>% arrange(group_close_to_gwas)
median_pc_alzheimer$group_close_to_gwas <- 1:25
median_pc_alzheimer <- median_pc_alzheimer %>% arrange(weighted_pc_groups)
pgs_df_temp <- pgs_df %>%
  filter(phenotype == "Alzheimer") %>%
  left_join(median_pc_alzheimer[, c(1, 2, 4)], by = "weighted_pc_groups")

median_pc <- get_median_pc_disease_2(non_pgs_df)
median_pc$group_close_to_gwas <- abs(median_pc$median_pc - 1)
median_pc <- median_pc %>% arrange(group_close_to_gwas)
median_pc$group_close_to_gwas <- 1:250
median_pc <- median_pc %>% arrange(weighted_pc_groups)
pgs_df_temp_2 <- pgs_df %>%
  filter(phenotype != "Alzheimer") %>%
  left_join(median_pc[, c(1, 2, 4)], by = "weighted_pc_groups")

pgs_df <- pgs_df_temp %>% rbind.data.frame(pgs_df_temp_2)
rm(pgs_df_temp, pgs_df_temp_2 )

pgs_df <- pgs_df %>%
  arrange(phenotype, weighted_pc_groups)

non_pgs_df %>% write_tsv("data/pgs_pred/group_non_pgs_df_disease.tsv")
pgs_df %>% write_tsv("data/pgs_pred/group_pgs_df_disease.tsv")
