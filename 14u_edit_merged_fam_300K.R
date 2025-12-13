library(tidyverse)
library(stringr)

fam_file <- read_delim('/scratch/06568/joyce_w/pgs_portability_questions/data/ukb_merged/ukb_merged.fam',
                       col_names = c('#FID', 'IID', 'V1', 'V2', 'V3', 'V4'),
                       delim = ' ', trim_ws = T)

nwb <- read_delim('data/ukb_populations/nwb_all_id_300K.txt', delim = ' ', trim_ws = T)

wb_gwas <- read_delim('data/ukb_populations/wb_gwas_id_300K.txt', delim = ' ', trim_ws = T)

wb_pred <- read_delim('data/ukb_populations/wb_pred_id_300K.txt', delim = ' ', trim_ws = T)

combined_labels <- bind_rows(
  wb_gwas %>% mutate(pop = 'GWAS'),
  wb_pred %>% mutate(pop = str_glue('pred{IID}') %>% as.character),
  nwb %>% mutate(pop = str_glue('pred{IID}') %>% as.character)
)

combined <- fam_file %>%
  left_join(combined_labels, by = c('#FID', 'IID'))

combined %>%
  select(-'#FID') %>%
  select('#FID' = pop, 'IID', starts_with('V')) %>%
  write_tsv('/scratch/06568/joyce_w/pgs_portability_questions/data/ukb_merged/ukb_merged_edited_300K.fam', col_names = F)
