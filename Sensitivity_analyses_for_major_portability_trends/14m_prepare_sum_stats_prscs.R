library(dplyr)
library(data.table)
library(readr)

pheno <- c("BMI", "Lymphocyte", "Height", "Eosinophil", "MCH", "MCV", "Monocyte", "Platelet", "RBC", 
           "WBC", "LDL", "Weight", "Triglycerides", "Cystatin_C", "Body_Fat_Perc")

for(i in pheno){
  for(chr in 1:22){
    file <- read_tsv(paste0("data/gwas_results/", i, ".chr", chr, i, ".glm.linear"))
    
    # Create a A2 columns based on A1
    file$A2 <- ifelse(file$ALT == file$A1, file$REF, file$ALT)
    
    colnames(file) <- c("CHR", "BP", "SNP", "REF", "ALT", "A1", "TEST", "OBS_CT", "BETA",
                        "BETA_SE", "T_STAT", "P", "ERRCODE", "A2")
    
    file <- file %>% select(SNP, A1, A2, BETA, P)
    
    file %>% write_tsv(paste0("data/gwas_results/", i, ".chr", chr, "_prscs.", i, ".glm.linear"))
  }
}
