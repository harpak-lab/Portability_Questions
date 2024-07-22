library(tidyverse)

# Import covariates
sex_df <- read_table("data/extracted_data_fields/genetic_sex.txt", col_names = T)
colnames(sex_df) = c("IID", "is_male")

age_df <- read_table("data/extracted_data_fields/age.txt", col_names = T)
colnames(age_df) = c("IID", "age")

pc_df <- read_table("data/extracted_data_fields/pc.txt", col_names=T)
cn = c()
for(i in 1:40){
  cn[i] = paste0("PC", i)
}

colnames(pc_df) = c("IID", cn)

df <- sex_df %>% full_join(age_df, by="IID") %>% full_join(pc_df, by="IID")
df <- df %>%
  filter(!is.na(age)) %>%
  rename(sex = is_male) %>%
  mutate(age_sq = age^2) %>%
  mutate(age_sex = age*sex) %>%
  mutate(age_sq_sex = age_sq * sex) %>%
  select(IID,sex, age, age_sq, age_sex, age_sq_sex, contains("PC"))

# Import population files
nwb_df <- read_delim("data/ukb_populations/nwb_all_id.txt", delim = " ", trim_ws = T)

wb_gwas_df <- read_delim("data/ukb_populations/wb_gwas_id.txt", delim = " ", trim_ws = T)
wb_pred_df <- read_delim("data/ukb_populations/wb_pred_id.txt", delim = " ", trim_ws = T)

combined_labels <- bind_rows(
  wb_gwas_df %>% mutate(pop = 'WB_GWAS'),
  wb_pred_df %>% mutate(pop = 'WB_pred'),
  nwb_df %>% mutate(pop = 'NWB')
)

df_pop <- df %>%
  left_join(combined_labels, by = c('IID'))

df_pop <- cbind.data.frame(`#FID` = df_pop$`#FID`, df_pop %>% select(-`#FID`))

# Export file
df_pop %>% write_tsv('data/ukb_merged/covar.tsv')
