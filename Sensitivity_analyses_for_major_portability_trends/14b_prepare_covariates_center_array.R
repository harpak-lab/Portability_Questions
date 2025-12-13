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

array_df <- read_table('data/extracted_data_fields/array_type.txt', col_names=T)

array_df$UKBiLEVEAX_b11 = ifelse(array_df$`22000-0.0` == -11, 1, 0)
array_df$UKBiLEVEAX_b10 = ifelse(array_df$`22000-0.0` == -10, 1, 0)
array_df$UKBiLEVEAX_b9 = ifelse(array_df$`22000-0.0` == -9, 1, 0)
array_df$UKBiLEVEAX_b8 = ifelse(array_df$`22000-0.0` == -8, 1, 0)
array_df$UKBiLEVEAX_b7 = ifelse(array_df$`22000-0.0` == -7, 1, 0)
array_df$UKBiLEVEAX_b6 = ifelse(array_df$`22000-0.0` == -6, 1, 0)
array_df$UKBiLEVEAX_b5 = ifelse(array_df$`22000-0.0` == -5, 1, 0)
array_df$UKBiLEVEAX_b4 = ifelse(array_df$`22000-0.0` == -4, 1, 0)
array_df$UKBiLEVEAX_b3 = ifelse(array_df$`22000-0.0` == -3, 1, 0)
array_df$UKBiLEVEAX_b2 = ifelse(array_df$`22000-0.0` == -2, 1, 0)
array_df$UKBiLEVEAX_b1 = ifelse(array_df$`22000-0.0` == -1, 1, 0)
array_df$Batch_b001 = ifelse(array_df$`22000-0.0` == 1, 1, 0)
array_df$Batch_b002 = ifelse(array_df$`22000-0.0` == 2, 1, 0)
array_df$Batch_b003 = ifelse(array_df$`22000-0.0` == 3, 1, 0)
array_df$Batch_b004 = ifelse(array_df$`22000-0.0` == 4, 1, 0)
array_df$Batch_b005 = ifelse(array_df$`22000-0.0` == 5, 1, 0)
array_df$Batch_b006 = ifelse(array_df$`22000-0.0` == 6, 1, 0)
array_df$Batch_b007 = ifelse(array_df$`22000-0.0` == 7, 1, 0)
array_df$Batch_b008 = ifelse(array_df$`22000-0.0` == 8, 1, 0)
array_df$Batch_b009 = ifelse(array_df$`22000-0.0` == 9, 1, 0)
array_df$Batch_b010 = ifelse(array_df$`22000-0.0` == 10, 1, 0)
array_df$Batch_b011 = ifelse(array_df$`22000-0.0` == 11, 1, 0)
array_df$Batch_b012 = ifelse(array_df$`22000-0.0` == 12, 1, 0)
array_df$Batch_b013 = ifelse(array_df$`22000-0.0` == 13, 1, 0)
array_df$Batch_b014 = ifelse(array_df$`22000-0.0` == 14, 1, 0)
array_df$Batch_b015 = ifelse(array_df$`22000-0.0` == 15, 1, 0)
array_df$Batch_b016 = ifelse(array_df$`22000-0.0` == 16, 1, 0)
array_df$Batch_b017 = ifelse(array_df$`22000-0.0` == 17, 1, 0)
array_df$Batch_b018 = ifelse(array_df$`22000-0.0` == 18, 1, 0)
array_df$Batch_b019 = ifelse(array_df$`22000-0.0` == 19, 1, 0)
array_df$Batch_b020 = ifelse(array_df$`22000-0.0` == 20, 1, 0)
array_df$Batch_b021 = ifelse(array_df$`22000-0.0` == 21, 1, 0)
array_df$Batch_b022 = ifelse(array_df$`22000-0.0` == 22, 1, 0)
array_df$Batch_b023 = ifelse(array_df$`22000-0.0` == 23, 1, 0)
array_df$Batch_b024 = ifelse(array_df$`22000-0.0` == 24, 1, 0)
array_df$Batch_b025 = ifelse(array_df$`22000-0.0` == 25, 1, 0)
array_df$Batch_b026 = ifelse(array_df$`22000-0.0` == 26, 1, 0)
array_df$Batch_b027 = ifelse(array_df$`22000-0.0` == 27, 1, 0)
array_df$Batch_b028 = ifelse(array_df$`22000-0.0` == 28, 1, 0)
array_df$Batch_b029 = ifelse(array_df$`22000-0.0` == 29, 1, 0)
array_df$Batch_b030 = ifelse(array_df$`22000-0.0` == 30, 1, 0)
array_df$Batch_b031 = ifelse(array_df$`22000-0.0` == 31, 1, 0)
array_df$Batch_b032 = ifelse(array_df$`22000-0.0` == 32, 1, 0)
array_df$Batch_b033 = ifelse(array_df$`22000-0.0` == 33, 1, 0)
array_df$Batch_b034 = ifelse(array_df$`22000-0.0` == 34, 1, 0)
array_df$Batch_b035 = ifelse(array_df$`22000-0.0` == 35, 1, 0)
array_df$Batch_b036 = ifelse(array_df$`22000-0.0` == 36, 1, 0)
array_df$Batch_b037 = ifelse(array_df$`22000-0.0` == 37, 1, 0)
array_df$Batch_b038 = ifelse(array_df$`22000-0.0` == 38, 1, 0)
array_df$Batch_b039 = ifelse(array_df$`22000-0.0` == 39, 1, 0)
array_df$Batch_b040 = ifelse(array_df$`22000-0.0` == 40, 1, 0)
array_df$Batch_b041 = ifelse(array_df$`22000-0.0` == 41, 1, 0)
array_df$Batch_b042 = ifelse(array_df$`22000-0.0` == 42, 1, 0)
array_df$Batch_b043 = ifelse(array_df$`22000-0.0` == 43, 1, 0)
array_df$Batch_b044 = ifelse(array_df$`22000-0.0` == 44, 1, 0)
array_df$Batch_b045 = ifelse(array_df$`22000-0.0` == 45, 1, 0)
array_df$Batch_b046 = ifelse(array_df$`22000-0.0` == 46, 1, 0)
array_df$Batch_b047 = ifelse(array_df$`22000-0.0` == 47, 1, 0)
array_df$Batch_b048 = ifelse(array_df$`22000-0.0` == 48, 1, 0)
array_df$Batch_b049 = ifelse(array_df$`22000-0.0` == 49, 1, 0)
array_df$Batch_b050 = ifelse(array_df$`22000-0.0` == 50, 1, 0)
array_df$Batch_b051 = ifelse(array_df$`22000-0.0` == 51, 1, 0)
array_df$Batch_b052 = ifelse(array_df$`22000-0.0` == 52, 1, 0)
array_df$Batch_b053 = ifelse(array_df$`22000-0.0` == 53, 1, 0)
array_df$Batch_b054 = ifelse(array_df$`22000-0.0` == 54, 1, 0)
array_df$Batch_b055 = ifelse(array_df$`22000-0.0` == 55, 1, 0)
array_df$Batch_b056 = ifelse(array_df$`22000-0.0` == 56, 1, 0)
array_df$Batch_b057 = ifelse(array_df$`22000-0.0` == 57, 1, 0)
array_df$Batch_b058 = ifelse(array_df$`22000-0.0` == 58, 1, 0)
array_df$Batch_b059 = ifelse(array_df$`22000-0.0` == 59, 1, 0)
array_df$Batch_b060 = ifelse(array_df$`22000-0.0` == 60, 1, 0)
array_df$Batch_b061 = ifelse(array_df$`22000-0.0` == 61, 1, 0)
array_df$Batch_b062 = ifelse(array_df$`22000-0.0` == 62, 1, 0)
array_df$Batch_b063 = ifelse(array_df$`22000-0.0` == 63, 1, 0)
array_df$Batch_b064 = ifelse(array_df$`22000-0.0` == 64, 1, 0)
array_df$Batch_b065 = ifelse(array_df$`22000-0.0` == 65, 1, 0)
array_df$Batch_b066 = ifelse(array_df$`22000-0.0` == 66, 1, 0)
array_df$Batch_b067 = ifelse(array_df$`22000-0.0` == 67, 1, 0)
array_df$Batch_b068 = ifelse(array_df$`22000-0.0` == 68, 1, 0)
array_df$Batch_b069 = ifelse(array_df$`22000-0.0` == 69, 1, 0)
array_df$Batch_b070 = ifelse(array_df$`22000-0.0` == 70, 1, 0)
array_df$Batch_b071 = ifelse(array_df$`22000-0.0` == 71, 1, 0)
array_df$Batch_b072 = ifelse(array_df$`22000-0.0` == 72, 1, 0)
array_df$Batch_b073 = ifelse(array_df$`22000-0.0` == 73, 1, 0)
array_df$Batch_b074 = ifelse(array_df$`22000-0.0` == 74, 1, 0)
array_df$Batch_b075 = ifelse(array_df$`22000-0.0` == 75, 1, 0)
array_df$Batch_b076 = ifelse(array_df$`22000-0.0` == 76, 1, 0)
array_df$Batch_b077 = ifelse(array_df$`22000-0.0` == 77, 1, 0)
array_df$Batch_b078 = ifelse(array_df$`22000-0.0` == 78, 1, 0)
array_df$Batch_b079 = ifelse(array_df$`22000-0.0` == 79, 1, 0)
array_df$Batch_b080 = ifelse(array_df$`22000-0.0` == 80, 1, 0)
array_df$Batch_b081 = ifelse(array_df$`22000-0.0` == 81, 1, 0)
array_df$Batch_b082 = ifelse(array_df$`22000-0.0` == 82, 1, 0)
array_df$Batch_b083 = ifelse(array_df$`22000-0.0` == 83, 1, 0)
array_df$Batch_b084 = ifelse(array_df$`22000-0.0` == 84, 1, 0)
array_df$Batch_b085 = ifelse(array_df$`22000-0.0` == 85, 1, 0)
array_df$Batch_b086 = ifelse(array_df$`22000-0.0` == 86, 1, 0)
array_df$Batch_b087 = ifelse(array_df$`22000-0.0` == 87, 1, 0)
array_df$Batch_b088 = ifelse(array_df$`22000-0.0` == 88, 1, 0)
array_df$Batch_b089 = ifelse(array_df$`22000-0.0` == 89, 1, 0)
array_df$Batch_b090 = ifelse(array_df$`22000-0.0` == 90, 1, 0)
array_df$Batch_b091 = ifelse(array_df$`22000-0.0` == 91, 1, 0)
array_df$Batch_b092 = ifelse(array_df$`22000-0.0` == 92, 1, 0)
array_df$Batch_b093 = ifelse(array_df$`22000-0.0` == 93, 1, 0)
array_df$Batch_b094 = ifelse(array_df$`22000-0.0` == 94, 1, 0)
array_df$Batch_b095 = ifelse(array_df$`22000-0.0` == 95, 1, 0)

colnames(array_df)[1] = "IID"

center_df <- read_table('data/extracted_data_fields/assessment_center.txt', col_names=T)
center_df <- center_df[, c(1:2)]

center_df$Barts = ifelse(center_df$`54-0.0` == 11012, 1, 0)
center_df$Birmingham = ifelse(center_df$`54-0.0` == 11021, 1, 0)
center_df$Bristol = ifelse(center_df$`54-0.0` == 11011, 1, 0)
center_df$Bury = ifelse(center_df$`54-0.0` == 11008, 1, 0)
center_df$Cardiff = ifelse(center_df$`54-0.0` == 11003, 1, 0)
center_df$Croydon = ifelse(center_df$`54-0.0` == 11020, 1, 0)
center_df$Edinburgh = ifelse(center_df$`54-0.0` == 11005, 1, 0)
center_df$Glasgow = ifelse(center_df$`54-0.0` == 11004, 1, 0)
center_df$Hounslow = ifelse(center_df$`54-0.0` == 11018, 1, 0)
center_df$Leeds = ifelse(center_df$`54-0.0` == 11010, 1, 0)
center_df$Liverpool = ifelse(center_df$`54-0.0` == 11016, 1, 0)
center_df$Manchester = ifelse(center_df$`54-0.0` == 11001, 1, 0)
center_df$Middlesborough = ifelse(center_df$`54-0.0` == 11017, 1, 0)
center_df$Newcastle = ifelse(center_df$`54-0.0` == 11009, 1, 0)
center_df$Nottingham = ifelse(center_df$`54-0.0` == 11013, 1, 0)
center_df$Oxford = ifelse(center_df$`54-0.0` == 11002, 1, 0)
center_df$Reading = ifelse(center_df$`54-0.0` == 11007, 1, 0)
center_df$Sheffield = ifelse(center_df$`54-0.0` == 11014, 1, 0)
center_df$Stockport = ifelse(center_df$`54-0.0` == 10003, 1, 0)
center_df$Stoke = ifelse(center_df$`54-0.0` == 11006, 1, 0)
center_df$Swansea = ifelse(center_df$`54-0.0` == 11022, 1, 0)
center_df$Wrexham = ifelse(center_df$`54-0.0` == 11023, 1, 0)

colnames(center_df)[1] = "IID"

df <- sex_df %>% full_join(age_df,by="IID") %>% full_join(pc_df,by="IID") %>% 
  full_join(array_df,by="IID") %>% full_join(center_df,by="IID")
df <- df %>%
  filter(!is.na(age)) %>%
  # Code sex as 0 = missing, 1 = female, 2 = male, as in plink .sample files
  rename(sex_covar = is_male) %>%
  mutate(age_sq = age^2) %>%
  mutate(age_sex = age*sex_covar) %>%
  mutate(age_sq_sex = age_sq * sex_covar) %>%
  select(IID, sex_covar, age, age_sq, age_sex, age_sq_sex, contains("PC"), 
         UKBiLEVEAX_b11:Batch_b094, Barts:Swansea)


# Import population files
nwb_df <- read_delim("data/ukb_populations/nwb_all_id.txt", delim = " ", trim_ws = T)

wb_gwas_df <- read_delim("data/ukb_populations/wb_gwas_id.txt", delim = " ", trim_ws = T)

combined_labels <- bind_rows(
  wb_gwas_df %>% mutate(pop = 'WB_GWAS'),
  nwb_df %>% mutate(pop = 'NWB')
)

df_pop <- df %>%
  left_join(combined_labels, by = c('IID'))

df_pop <- cbind.data.frame(`#FID` = df_pop$`#FID`, df_pop %>% select(-`#FID`))

# Export file
df_pop %>% write_tsv('data/ukb_merged/covar_array_center.tsv')
