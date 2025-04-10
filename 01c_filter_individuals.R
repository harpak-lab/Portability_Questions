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
not_related = read_table("data/extracted_data_fields/used_in_pca.txt")

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

keep_id = cbind.data.frame(`#FID` = filtered_4, IID = filtered_4)
keep_id %>% write.table("data/ukb_populations/keep_id.txt", 
                        row.names = F, col.names = T, quote = F, sep = " ")

# Data fields for identifying WB
ethnic_background <- read_table("data/extracted_data_fields/ethnic_background.txt")

wb <- filtered_4[filtered_4 %in% 
                   ethnic_background$eid[ethnic_background$`21000-0.0` == 1001]]
nwb <- filtered_4[filtered_4 %notin% wb]

wb <- cbind.data.frame(`#FID` = wb, IID = wb)
nwb <- cbind.data.frame(`#FID` = nwb, IID = nwb)

wb %>% write.table("data/ukb_populations/wb_all_id.txt", 
                   row.names = F, col.names = T, quote = F, sep = " ")

nwb %>% write.table("data/ukb_populations/nwb_all_id.txt", 
                    row.names = F, col.names = T, quote = F, sep = " ")

# For GWAS, first remove relatives from WB
gwas <- wb %>% filter(IID %in% not_related$eid[!is.na(not_related$`22020-0.0`) & not_related$`22020-0.0` == 1])

# Sample 350K WB for GWAS
set.seed(3)
gwas <- sample_n(gwas, 350000)

gwas <- gwas %>% arrange(IID)
gwas %>% write.table("data/ukb_populations/wb_gwas_id.txt",
                     sep = " ", col.names = T, row.names = F, quote = F)

wb_pred <- wb %>% filter(IID %notin% gwas$IID)
wb_pred %>% write.table("data/ukb_populations/wb_pred_id.txt",
                        sep = " ", col.names = T, row.names = F, quote = F)