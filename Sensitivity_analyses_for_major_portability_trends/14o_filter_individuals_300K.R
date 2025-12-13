library(dplyr)
library(tidyr)
library(readr)

`%notin%` <- Negate(`%in%`)

# Read in ID files
keep_id <- read_delim("data/ukb_populations/keep_id.txt", delim = " ", trim_ws = T)
wb <- read_delim("data/ukb_populations/wb_gwas_id.txt", delim = " ", trim_ws = T)
nwb <- read_delim("data/ukb_populations/nwb_all_id.txt", delim = " ", trim_ws = T)

# Randomly select 300K WB for GWAS
set.seed(3)
gwas <- sample(wb, size = 300000, replace = F)
gwas <- cbind.data.frame("#FID" = gwas, IID = gwas)
gwas <- gwas %>% arrange(IID)

# Use the rest of WB and all NWB for prediction
prediction_wb <- wb[wb %notin% gwas$IID]
prediction_wb <- cbind.data.frame("#FID" = prediction_wb, "IID" = prediction_wb)
prediction_nwb <- cbind.data.frame("#FID" = nwb, "IID" = nwb)

gwas %>% write.table("data/ukb_populations/wb_gwas_id_300K.txt", 
                   row.names = F, col.names = T, quote = F, sep = " ")

prediction_wb %>% write.table("data/ukb_populations/wb_pred_id_300K.txt", 
                   row.names = F, col.names = T, quote = F, sep = " ")

nwb %>% write.table("data/ukb_populations/nwb_all_id_300K.txt", 
                    row.names = F, col.names = T, quote = F, sep = " ")
