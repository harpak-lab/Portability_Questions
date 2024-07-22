library(tidyr)
library(dplyr)
library(readr)

# Prepare covariates file
covar <- read_tsv("data/ukb_merged/covar.tsv")
covar <- covar %>% select(c(`#FID`, IID, sex, age, age_sq, age_sex, age_sq_sex))

close_pc <- read_tsv("data/pca/close_pca.eigenvec")
far_pc <- read_tsv("data/pca/far_pca.eigenvec")

covar_close <- covar %>% right_join(close_pc[, -1], by = "IID")
covar_far <- covar %>% right_join(far_pc[, -1], by = "IID")

covar_close %>% write_tsv("data/ukb_merged/covar_close.tsv")
covar_far %>% write_tsv("data/ukb_merged/covar_far.tsv")