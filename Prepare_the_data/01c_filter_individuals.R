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
not_related <- read_table("data/extracted_data_fields/used_in_pca.txt")
caucasion <- read_table("data/extracted_data_fields/caucasian.txt")
withdrew <- read_table("data/ukb_meta/w61666_20241217.csv", col_names = F)

unfiltered <- unfiltered[unfiltered %notin% withdrew$X1]

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

# Remove relatives
filtered_5 <- filtered_4[filtered_4 %in% not_related$eid[not_related$`22020-0.0` == 1]]


keep_id = cbind.data.frame(`#FID` = filtered_4, IID = filtered_4)
keep_id %>% write.table("data/ukb_populations/keep_id.txt", 
                        row.names = F, col.names = T, quote = F, sep = " ")

# Separate WB and NWB
wb <- filtered_5[filtered_5 %in% caucasion$eid[caucasion$`22006-0.0` == 1]]
nwb <- filtered_5[filtered_5 %notin% caucasion$eid[caucasion$`22006-0.0` == 1]]

wb <- cbind.data.frame(`#FID` = wb, IID = wb)
nwb <- cbind.data.frame(`#FID` = nwb, IID = nwb)

# Use all WB for GWAS
wb %>% write.table("data/ukb_populations/wb_gwas_id.txt", 
                   row.names = F, col.names = T, quote = F, sep = " ")

# Use all NWB for prediction
nwb %>% write.table("data/ukb_populations/nwb_all_id.txt", 
                    row.names = F, col.names = T, quote = F, sep = " ")
