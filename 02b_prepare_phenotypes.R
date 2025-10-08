library(tidyverse)
library(tibble)

# Table of each trait's name, field ID, and data type
traits_info <- cbind.data.frame(
  trait = c("Height", "Cystatin_C", "Platelet", "MCV", "Weight", 
            "MCH", "BMI", "RBC", "Body_Fat_Perc", "Monocyte", 
            "Triglycerides", "Lymphocyte", "WBC", "Eosinophil", "LDL"),
  ukb_field = c("50-0.0", "30720-0.0", "30080-0.0", "30040-0.0", "21002-0.0",
                "30050-0.0", "21001-0.0", "30010-0.0", "23099-0.0", "30130-0.0",
                "30870-0.0", "30120-0.0", "30000-0.0", "30150-0.0", "30780-0.0"),
  dtype = "d") %>%
  add_row(trait = "eid", ukb_field = "eid", dtype = "c")

# Load only the 15 columns of interest
field_to_dtype <- traits_info %>%
  select(ukb_field, dtype) %>%
  deframe %>%
  as.list

raw_phenotypes_df <- read_csv('/work/06568/joyce_w/stampede2/software/ukbconv/ukb45020.csv',
                              col_types = do.call(cols_only, field_to_dtype))

# Import population files
nwb_df <- read_delim("data/ukb_populations/nwb_all_id.txt", delim = " ", trim_ws = T)

wb_gwas_df <- read_delim("data/ukb_populations/wb_gwas_id.txt", delim = " ", trim_ws = T)

combined_labels <- bind_rows(
  wb_gwas_df %>% mutate(pop = 'WB_GWAS'),
  nwb_df %>% mutate(pop = 'NWB')
)
combined_labels$IID <- as.character(combined_labels$IID)
combined_labels$`#FID` <- as.character(combined_labels$`#FID`)

pheno_df <- raw_phenotypes_df %>%
  # Rename fields to trait names
  pivot_longer(-eid, names_to = 'ukb_field') %>%
  inner_join(traits_info %>% select(-dtype), by = 'ukb_field') %>%
  pivot_wider(id_cols = eid, names_from = trait, values_from = value) %>%
  rename(IID = eid) %>%
  left_join(combined_labels,by=c('IID'))

pheno_df %>% select('#FID', IID, all_of(traits_info$trait[1:15])) %>%
  write_tsv('data/phenotypes/phenotypes.tsv')
