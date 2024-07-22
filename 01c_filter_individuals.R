library(dplyr)
library(tidyr)
library(readr)

`%notin%` <- Negate(`%in%`)

unfiltered <- read_table("data/extracted_data_fields/sex_chrom_aneuploidy.txt")
unfiltered <- unfiltered$eid

# Data fields for filtering
sex_chrom_aneuploidy <- read_table("data/extracted_data_fields/sex_chrom_aneuploidy.txt")
reported_sex <- read_table("data/extracted_data_fields/reported_sex.txt")
gen_sex <- read_table("data/extracted_data_fields/genetic_sex.txt")
outlier_missingness <- read_table("data/extracted_data_fields/outlier_heterozygosity.txt")
missingness <- read_table("data/extracted_data_fields/genotype_missingness.txt")
relatives = read_table("data/extracted_data_fields/genetic_kinship.txt")

# Remove those with sex chromosome aneuploidy
filtered_1 <- unfiltered[unfiltered %notin% 
                           sex_chrom_aneuploidy$eid[sex_chrom_aneuploidy$`22019-0.0` == 1]]

sex_diff = reported_sex %>% left_join(gen_sex, by = "eid")
# Keep those whose self-reported sex is the same as sex determined from genotyping analysis
filtered_2 <- filtered_1[filtered_1 %in% 
                           sex_diff$eid[sex_diff$`31-0.0` == sex_diff$`22001-0.0`]]

# Remove those who are outliers in heterozygosity or genotype missingness
filtered_3 <- filtered_2[filtered_2 %notin% 
                          outlier_missingness$eid[outlier_missingness$`22027-0.0` == 1]]

# Keep those with a genotype missingness less than or equal to 2%
filtered_4 <- filtered_3[filtered_3 %in% 
                           missingness$eid[missingness$`22005-0.0` <= 0.02]]

# Keep those who have no kinship found or has at least one kinship identified, 
# but not those with 10 or more 3rd degree relatives identified
filtered_5 <- filtered_4[filtered_4 %in% 
                           relatives$eid[relatives$`22021-0.0` == 0 | relatives$`22021-0.0` == 1]]

keep_id = cbind.data.frame(`#FID` = filtered_5, IID = filtered_5)
keep_id %>% write.table("data/ukb_populations/keep_id.txt", 
                        row.names = F, col.names = T, quote = F, sep = " ")

# Data fields for identifying WB
ethnic_background <- read_table("data/extracted_data_fields/ethnic_background.txt")
caucasian <- read_table("data/extracted_data_fields/caucasian.txt")

wb <- filtered_5[filtered_5 %in% 
                   ethnic_background$eid[ethnic_background$`21000-0.0` == 1001]]
wb <- filtered_5[filtered_5 %in% 
                   caucasian$eid[caucasian$`22006-0.0` == 1]]
nwb <- filtered_5[filtered_5 %notin% wb]

wb <- cbind.data.frame(`#FID` = wb, IID = wb)
nwb <- cbind.data.frame(`#FID` = nwb, IID = nwb)

wb %>% write.table("data/ukb_populations/wb_all_id.txt", 
                   row.names = F, col.names = T, quote = F, sep = " ")

nwb %>% write.table("data/ukb_populations/nwb_all_id.txt", 
                    row.names = F, col.names = T, quote = F, sep = " ")

# Sample 350K WB for GWAS
set.seed(3)
wb_rand <- sample_n(wb, 350000)

wb_rand <- wb_rand %>% arrange(IID)
wb_rand %>% write.table("data/ukb_populations/wb_gwas_id.txt",
                         sep = " ", col.names = T, row.names = F, quote = F)

wb_pred <- wb %>% filter(IID %notin% wb_rand$IID)
wb_pred %>% write.table("data/ukb_populations/wb_pred_id.txt",
                        sep = " ", col.names = T, row.names = F, quote = F)