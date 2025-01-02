library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(stringr)

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
  covar_df <- covar_df %>% select(-c(8:47))
  
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
                        'data/ukb_populations/wb_gwas_id.txt', 
                        'data/ukb_populations/wb_pred_id.txt')
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
  
  for (file in list.files(path = 'data/pgs',
                          pattern = '[a-zA-Z]+_[0-9]_scores.sscore', full.names = T)) {
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
non_pgs_df <- load_non_pgs_df(num_bins = 500)
non_pgs_df <- non_pgs_df %>% filter(!is.na(weighted_pc_groups))
non_pgs_df <- non_pgs_df %>%
  mutate(array_type = as.factor(array_type))

# Data frame with PGS information
pgs_df  <- make_pgs_evaluation_df(non_pgs_df, "weighted_pc_groups")
pgs_df <- pgs_df %>% 
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

# Order the bins by how close to genetic distance = 1 (mean distance to the GWAS centroidof the GWAS group)
median_pc <- get_median_pc(non_pgs_df)
median_pc$group_close_to_gwas <- abs(median_pc$median_pc - 1)
median_pc <- median_pc %>% arrange(group_close_to_gwas)
median_pc$group_close_to_gwas <- 1:500
median_pc <- median_pc %>% arrange(weighted_pc_groups)
non_pgs_df <- non_pgs_df %>% left_join(median_pc[, c(1:2, 4)], by = "weighted_pc_groups")
pgs_df <- pgs_df %>% left_join(median_pc[, c(1:2, 4)], by = "weighted_pc_groups")

non_pgs_df %>% write_tsv("data/pgs_pred/group_non_pgs_df.tsv")
pgs_df %>% write_tsv("data/pgs_pred/group_pgs_df.tsv")