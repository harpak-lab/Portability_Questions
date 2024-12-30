library(dplyr)
library(tidyr)
library(readr)
library(data.table)

gwas_id <- read_table('data/ukb_populations/wb_gwas_id.txt')
icd_main <- read_table('data/extracted_data_fields/icd_10_main.txt')
icd_secondary <- read_table('data/extracted_data_fields/icd_10_secondary.txt')

# Mark as 1 if someone is positive for the disease
icd_main$alzheimer_main <- apply(icd_main, 1, function(row) {
  # Check if all columns except the first are NA
  if (all(is.na(row[-1]))) {
    return(NA)  # Set as NA if all columns except the first are NA
  } else {
    # Check if "G30" appears in any non-NA column
    return(ifelse(any(grepl("^(G30|F00)", row[!is.na(row)])), 1, 0))
  }
})

icd_main$t2d_main <- apply(icd_main, 1, function(row) {
  # Check if all columns except the first and the last are NA
  if (all(is.na(row[-c(1, 68)]))) {
    return(NA)  # Set as NA if all columns except the first are NA
  } else {
    # Check if "E11" appears in any non-NA column
    return(ifelse(any(grepl("^E11", row[!is.na(row)])), 1, 0))
  }
})

icd_main$asthma_main <- apply(icd_main, 1, function(row) {
  # Check if all columns except the first and the last are NA
  if (all(is.na(row[-c(1, 69)]))) {
    return(NA)  # Set as NA if all columns except the first are NA
  } else {
    # Check if "J45" appears in any non-NA column
    return(ifelse(any(grepl("^J45", row[!is.na(row)])), 1, 0))
  }
})

icd_combined <- icd_main %>%
  select(eid, alzheimer_main, t2d_main, asthma_main)

# Mark as 1 if someone is positive for the disease
icd_secondary$alzheimer_sec <- apply(icd_secondary, 1, function(row) {
  # Check if all columns except the first are NA
  if (all(is.na(row[-1]))) {
    return(NA)  # Set as NA if all columns except the first are NA
  } else {
    # Check if "G30" appears in any non-NA column
    return(ifelse(any(grepl("^(G30|F00)", row[!is.na(row)])), 1, 0))
  }
})

icd_secondary$t2d_sec <- apply(icd_secondary, 1, function(row) {
  # Check if all columns except the first and the last are NA
  if (all(is.na(row[-c(1, 186)]))) {
    return(NA)  # Set as NA if all columns except the first are NA
  } else {
    # Check if "E11" appears in any non-NA column
    return(ifelse(any(grepl("^E11", row[!is.na(row)])), 1, 0))
  }
})

icd_secondary$asthma_sec <- apply(icd_secondary, 1, function(row) {
  # Check if all columns except the first and the last are NA
  if (all(is.na(row[-c(1, 187)]))) {
    return(NA)  # Set as NA if all columns except the first are NA
  } else {
    # Check if "J45" appears in any non-NA column
    return(ifelse(any(grepl("^J45", row[!is.na(row)])), 1, 0))
  }
})

icd_combined <- icd_combined %>%
  left_join(icd_secondary[c(1, 186:188)], by = "eid")

icd_combined$alzheimer <- mapply(function(x, y) {
  if (is.na(x) & is.na(y)) {
    return(NA)              # Keep as NA if both are NA
  } else if (is.na(x)) {
    return(y)               # Use y if x is NA
  } else if (is.na(y)) {
    return(x)               # Use x if y is NA
  } else {
    return(min(x + y, 1))   # Cap the sum at 1 if both are non-NA
  }
}, icd_combined$alzheimer_main, icd_combined$alzheimer_sec)

icd_combined$t2d <- mapply(function(x, y) {
  if (is.na(x) & is.na(y)) {
    return(NA)              # Keep as NA if both are NA
  } else if (is.na(x)) {
    return(y)               # Use y if x is NA
  } else if (is.na(y)) {
    return(x)               # Use x if y is NA
  } else {
    return(min(x + y, 1))   # Cap the sum at 1 if both are non-NA
  }
}, icd_combined$t2d_main, icd_combined$t2d_sec)

icd_combined$asthma <- mapply(function(x, y) {
  if (is.na(x) & is.na(y)) {
    return(NA)              # Keep as NA if both are NA
  } else if (is.na(x)) {
    return(y)               # Use y if x is NA
  } else if (is.na(y)) {
    return(x)               # Use x if y is NA
  } else {
    return(min(x + y, 1))   # Cap the sum at 1 if both are non-NA
  }
}, icd_combined$asthma_main, icd_combined$asthma_sec)

icd_combined <- icd_combined %>% select(eid, alzheimer, t2d, asthma)
colnames(icd_combined) <- c("IID", "Alzheimer", "T2D", "Asthma")
icd_combined <- cbind.data.frame(`#FID` = icd_combined$IID, icd_combined)

icd_combined %>% write_tsv('data/phenotypes/phenotypes_disease.tsv')
