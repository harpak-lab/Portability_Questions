library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)
library(scales)
library(grid)
library(gridExtra)
library(stringr)
library(ggallin)

# Phenotypes
pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
           "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
           "Body_Fat_Perc")

# Read in gindividual PGS file
ind_pgs_df <- read_tsv("data/pgs_pred/ind_pgs_df.tsv") %>%
  select(-(16:34))

# Read in the predictors
townsend_index <- read_tsv("data/extracted_data_fields/townsend.txt")
townsend_index <- townsend_index %>%
  select(`189-0.0`) %>%
  separate(col = `189-0.0`, into = c("IID", "townsend"), sep = "\t")
townsend_index$IID <- as.numeric(townsend_index$IID)
townsend_index$townsend <- as.numeric(townsend_index$townsend)

ind_pgs_df <- ind_pgs_df %>%
  left_join(townsend_index, by = c("IID"))

income <- read.table("data/extracted_data_fields/income.txt", sep = "_", header = T)
income <- income %>%
  separate(col = 1, into = c("X", "IID", paste0("income_", 1:5)), sep = "\t") %>%
  select(-X) %>%
  select(IID, income_1) %>%
  mutate(IID = as.numeric(IID),
         income_1 = as.numeric(income_1))

income$income <- ifelse(income$income_1 == 1, "less_18", 
                       ifelse(income$income_1 == 2, "range_18_31",
                              ifelse(income$income_1 == 3, "range_31_52",
                                     ifelse(income$income_1 == 4, "range_52_100", "more_100"))))

income$income <- factor(income$income, levels = c("less_18", "range_18_31", "range_31_52",
                                                 "range_52_100", "more_100"))
income <- income %>% select(-income_1)

ind_pgs_df <- ind_pgs_df %>%
  left_join(income, by = c("IID"))

temp <- ind_pgs_df %>%
  filter(phenotype == "Height")
knots_pc_dist <- c(min(temp$pc_dist) + (max(temp$pc_dist) - min(temp$pc_dist)) / 5, 
                   min(temp$pc_dist) + (max(temp$pc_dist) - min(temp$pc_dist)) / 5 * 2, 
                   min(temp$pc_dist) + (max(temp$pc_dist) - min(temp$pc_dist)) / 5 * 3, 
                   min(temp$pc_dist) + (max(temp$pc_dist) - min(temp$pc_dist)) / 5 * 4,
                   max(temp$pc_dist))

temp <- temp %>%
  na.omit(townsend)
knots_townsend <- c(min(temp$townsend) + (max(temp$townsend) - min(temp$townsend)) / 5, 
                    min(temp$townsend) + (max(temp$townsend) - min(temp$townsend)) / 5 * 2, 
                    min(temp$townsend) + (max(temp$townsend) - min(temp$townsend)) / 5 * 3, 
                    min(temp$townsend) + (max(temp$townsend) - min(temp$townsend)) / 5 * 4,
                    max(temp$townsend))

plot_pred_dist <- function(pgs_df, trait, knots){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    group_by(phenotype) %>%
    na.omit(pc_dist, pred_error)
  
  # Plot prediction accuracy relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    filter(group_close_to_gwas <= 50) %>%
    dplyr::summarise(ref_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_error) %>%
    select(-ref_error) %>%
    ungroup()
  
  # Get the means
  plot_df$tile <- ifelse(plot_df$pc_dist <= knots[1], 1,
                        ifelse(plot_df$pc_dist <= knots[2], 2,
                               ifelse(plot_df$pc_dist <= knots[3], 3,
                                      ifelse(plot_df$pc_dist <= knots[4], 4, 5))))
  plot_df$tile <- factor(plot_df$tile)
  mean_df <- plot_df %>%
    dplyr::group_by(tile) %>%
    dplyr::summarize(median_pc = median(pc_dist, na.rm = T),
                     mean_pred = mean(relative_performance),
                     count = n(),
                     sd = sd(relative_performance),
                     se = sd / sqrt(count))
  
  if(trait == "Height"){
    c <- "#E69F00"
    val <- c("#ffc034", "#ffb000", "#E69F00", "#b37c00", "#9a6a00")
  } else if(trait == "Weight"){
    c <- "#40B0A6"
    val <- c("#71cbc3", "#4cbeb4", "#40B0A6", "#328b83", "#2c7871")
  } else if(trait == "LDL") {
    c <- "#DC3220"
    val <- c("#e86e61", "#e5594a", "#DC3220", "#af281a", "#992316")
  } else if(trait == "WBC"){
    c <- "#0072B2"
    val <- c("#00a3ff", "#0082cc", "#0072b2", "#00517f", "#004166")
  } else if(trait == 'BMI'){
    c <- "#9e79cc"
    val <- c("#b99eda", "#ab8bd3", "#9e79cc", "#8354be", "#8354be")
  } else if(trait == 'Body_Fat_Perc'){
    c <- "#cc9e79"
    val <- c("#dab99e", "#d3ab8b", "#cc9e79", "#be8354", "#b47645")
  } else if(trait == 'Monocyte'){
    c <- "#CC79A7"
    val <- c("#e1b0cb", "#d38bb3", "#CC79A7", "#be548f", "#b44582")
  } else if(trait == 'Lymphocyte'){
    c <- "#009E73"
    val <- c("#00d198", "#00b886", "#009E73", "#006b4e", "#00523b")
  } else if(trait == 'Eosinophil'){
    c <- "#79ccc7"
    val <- c("#9edad6", "#8bd3cf", "#79ccc7", "#54beb8", "#45b4ad")
  } else if(trait == 'MCV'){
    c <- "#79cc9e"
    val <- c("#9edab9", "#8bd3ab", "#79cc9e", "#54be83", "#45b476")
  } else if(trait == 'MCH'){
    c <- "#9e002b"
    val <- c("#eb0040", "#d10039", "#9e002b", "#6b001d", "#520016")
  } else if(trait == 'RBC'){
    c <- "#b20072"
    val <- c("#ff00a3", "#e50093", "#b20072", "#7f0051", "#660041")
  } else if(trait == 'Cystatin_C'){
    c <- "#a7cc79"
    val <- c("#bfda9e", "#b3d38b", "#a7cc79", "#8fbe54", "#82b445")
  } else if(trait == 'Platelet'){
    c <- "#a640b0"
    val <- c("#c371cb", "#bb5ec5", "#a640b0", "#83328b", "#712c78")
  } else if(trait == 'Triglycerides'){
    c <- "#56B4E9"
    val <- c("#9ad2f2", "#83c8ef", "#56B4E9", "#29a0e3", "#1c93d7")
  }
  
  plot <- ggplot(mean_df, aes(x = median_pc, y = mean_pred)) +
    geom_errorbar(aes(ymin = mean_pred - se, ymax = mean_pred + se, color = tile),
                  width = 5, linewidth = 1) +
    geom_point(aes(color = tile), size = 3) +
    geom_line(size = 0.5, color = c) +
    scale_x_continuous(breaks = c(0, 40, 80, 120, 160), expand = c(0, 0)) +
    scale_color_manual(values = val) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica"),
          legend.position = "none") +
    scale_y_reverse()
  
  return(plot)
}

# Make the plots
plot_pred_dist_height <- plot_pred_dist(ind_pgs_df, "Height", knots_pc_dist)
plot_pred_dist_weight <- plot_pred_dist(ind_pgs_df, "Weight", knots_pc_dist)
plot_pred_dist_ldl <- plot_pred_dist(ind_pgs_df, "LDL", knots_pc_dist)
plot_pred_dist_wbc <- plot_pred_dist(ind_pgs_df, "WBC", knots_pc_dist)
plot_pred_dist_bmi <- plot_pred_dist(ind_pgs_df, "BMI", knots_pc_dist)
plot_pred_dist_body_fat_perc <- plot_pred_dist(ind_pgs_df, "Body_Fat_Perc", knots_pc_dist)
plot_pred_dist_monocyte <- plot_pred_dist(ind_pgs_df, "Monocyte", knots_pc_dist)
plot_pred_dist_lymphocyte <- plot_pred_dist(ind_pgs_df, "Lymphocyte", knots_pc_dist)
plot_pred_dist_eosinophil <- plot_pred_dist(ind_pgs_df, "Eosinophil", knots_pc_dist)
plot_pred_dist_mcv <- plot_pred_dist(ind_pgs_df, "MCV", knots_pc_dist)
plot_pred_dist_mch <- plot_pred_dist(ind_pgs_df, "MCH", knots_pc_dist)
plot_pred_dist_rbc <- plot_pred_dist(ind_pgs_df, "RBC", knots_pc_dist)
plot_pred_dist_cystatin_c <- plot_pred_dist(ind_pgs_df, "Cystatin_C", knots_pc_dist)
plot_pred_dist_platelet <- plot_pred_dist(ind_pgs_df, "Platelet", knots_pc_dist)
plot_pred_dist_triglycerides <- plot_pred_dist(ind_pgs_df, "Triglycerides", knots_pc_dist)

plot_pred_townsend <- function(pgs_df, trait, knots){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    na.omit(townsend, pred_error)
  
  # Plot prediction accuracy relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 50) %>%
    dplyr::summarise(ref_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_error) %>%
    select(-ref_error) %>%
    ungroup()
  
  # Get the means
  plot_df$tile <- ifelse(plot_df$townsend <= knots[1], 1,
                        ifelse(plot_df$townsend <= knots[2], 2,
                               ifelse(plot_df$townsend <= knots[3], 3,
                                      ifelse(plot_df$townsend <= knots[4], 4, 5))))
  plot_df$tile <- factor(plot_df$tile)
  mean_df <- plot_df %>%
    dplyr::group_by(tile) %>%
    dplyr::summarize(median_townsend = median(townsend, na.rm = T),
                     mean_pred = mean(relative_performance),
                     count = n(),
                     sd = sd(relative_performance),
                     se = sd / sqrt(count))
  
  if(trait == "Height"){
    c <- "#E69F00"
    val <- c("#ffc034", "#ffb000", "#E69F00", "#b37c00", "#9a6a00")
  } else if(trait == "Weight"){
    c <- "#40B0A6"
    val <- c("#71cbc3", "#4cbeb4", "#40B0A6", "#328b83", "#2c7871")
  } else if(trait == "LDL") {
    c <- "#DC3220"
    val <- c("#e86e61", "#e5594a", "#DC3220", "#af281a", "#992316")
  } else if(trait == "WBC"){
    c <- "#0072B2"
    val <- c("#00a3ff", "#0082cc", "#0072b2", "#00517f", "#004166")
  } else if(trait == 'BMI'){
    c <- "#9e79cc"
    val <- c("#b99eda", "#ab8bd3", "#9e79cc", "#8354be", "#8354be")
  } else if(trait == 'Body_Fat_Perc'){
    c <- "#cc9e79"
    val <- c("#dab99e", "#d3ab8b", "#cc9e79", "#be8354", "#b47645")
  } else if(trait == 'Monocyte'){
    c <- "#CC79A7"
    val <- c("#e1b0cb", "#d38bb3", "#CC79A7", "#be548f", "#b44582")
  } else if(trait == 'Lymphocyte'){
    c <- "#009E73"
    val <- c("#00d198", "#00b886", "#009E73", "#006b4e", "#00523b")
  } else if(trait == 'Eosinophil'){
    c <- "#79ccc7"
    val <- c("#9edad6", "#8bd3cf", "#79ccc7", "#54beb8", "#45b4ad")
  } else if(trait == 'MCV'){
    c <- "#79cc9e"
    val <- c("#9edab9", "#8bd3ab", "#79cc9e", "#54be83", "#45b476")
  } else if(trait == 'MCH'){
    c <- "#9e002b"
    val <- c("#eb0040", "#d10039", "#9e002b", "#6b001d", "#520016")
  } else if(trait == 'RBC'){
    c <- "#b20072"
    val <- c("#ff00a3", "#e50093", "#b20072", "#7f0051", "#660041")
  } else if(trait == 'Cystatin_C'){
    c <- "#a7cc79"
    val <- c("#bfda9e", "#b3d38b", "#a7cc79", "#8fbe54", "#82b445")
  } else if(trait == 'Platelet'){
    c <- "#a640b0"
    val <- c("#c371cb", "#bb5ec5", "#a640b0", "#83328b", "#712c78")
  } else if(trait == 'Triglycerides'){
    c <- "#56B4E9"
    val <- c("#9ad2f2", "#83c8ef", "#56B4E9", "#29a0e3", "#1c93d7")
  }
  
  plot <- ggplot(mean_df, aes(x = median_townsend, y = mean_pred)) +
    geom_errorbar(aes(ymin = mean_pred - se, ymax = mean_pred + se, color = tile),
                  width = 0.3, linewidth = 1) +
    geom_point(aes(color = tile), size = 3) +
    geom_line(size = 0.5, color = c) +
    scale_x_continuous(breaks = c(-2.5, 0.0, 2.5, 5.0, 7.5), 
                       labels = c("-2.5", "0.0", "2.5", "5.0", "7.5"), expand = c(0, 0)) +
    scale_color_manual(values = val) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica"),
          legend.position = "none") +
    scale_y_reverse()
  
  return(plot)
}

# Make the plots
plot_pred_townsend_height <- plot_pred_townsend(ind_pgs_df, "Height", knots_townsend)
plot_pred_townsend_weight <- plot_pred_townsend(ind_pgs_df, "Weight", knots_townsend)
plot_pred_townsend_ldl <- plot_pred_townsend(ind_pgs_df, "LDL", knots_townsend)
plot_pred_townsend_wbc <- plot_pred_townsend(ind_pgs_df, "WBC", knots_townsend)
plot_pred_townsend_bmi <- plot_pred_townsend(ind_pgs_df, "BMI", knots_townsend)
plot_pred_townsend_body_fat_perc <- plot_pred_townsend(ind_pgs_df, "Body_Fat_Perc", knots_townsend)
plot_pred_townsend_monocyte <- plot_pred_townsend(ind_pgs_df, "Monocyte", knots_townsend)
plot_pred_townsend_lymphocyte <- plot_pred_townsend(ind_pgs_df, "Lymphocyte", knots_townsend)
plot_pred_townsend_eosinophil <- plot_pred_townsend(ind_pgs_df, "Eosinophil", knots_townsend)
plot_pred_townsend_mcv <- plot_pred_townsend(ind_pgs_df, "MCV", knots_townsend)
plot_pred_townsend_mch <- plot_pred_townsend(ind_pgs_df, "MCH", knots_townsend)
plot_pred_townsend_rbc <- plot_pred_townsend(ind_pgs_df, "RBC", knots_townsend)
plot_pred_townsend_cystatin_c <- plot_pred_townsend(ind_pgs_df, "Cystatin_C", knots_townsend)
plot_pred_townsend_platelet <- plot_pred_townsend(ind_pgs_df, "Platelet", knots_townsend)
plot_pred_townsend_triglycerides <- plot_pred_townsend(ind_pgs_df, "Triglycerides", knots_townsend)

plot_pred_income <- function(pgs_df, trait){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    na.omit(income, pred_error)
  
  # Plot prediction accuracy relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 50) %>%
    dplyr::summarise(ref_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_error) %>%
    select(-ref_error) %>%
    ungroup()
  
  # Get the means
  plot_df$tile <- factor(plot_df$income, labels = c("<18K", "18K-31K", "31K-52K", 
                                                   "52K-100K", ">100K"))
  mean_df <- plot_df %>%
    dplyr::group_by(tile) %>%
    dplyr::summarize(mean_pred = mean(relative_performance),
                     count = n(),
                     sd = sd(relative_performance),
                     se = sd / sqrt(count))
  mean_df$phenotype <- trait
  
  if(trait == "Height"){
    c <- "#E69F00"
    val <- c("#ffc034", "#ffb000", "#E69F00", "#b37c00", "#9a6a00")
  } else if(trait == "Weight"){
    c <- "#40B0A6"
    val <- c("#71cbc3", "#4cbeb4", "#40B0A6", "#328b83", "#2c7871")
  } else if(trait == "LDL") {
    c <- "#DC3220"
    val <- c("#e86e61", "#e5594a", "#DC3220", "#af281a", "#992316")
  } else if(trait == "WBC"){
    c <- "#0072B2"
    val <- c("#00a3ff", "#0082cc", "#0072b2", "#00517f", "#004166")
  } else if(trait == 'BMI'){
    c <- "#9e79cc"
    val <- c("#b99eda", "#ab8bd3", "#9e79cc", "#8354be", "#8354be")
  } else if(trait == 'Body_Fat_Perc'){
    c <- "#cc9e79"
    val <- c("#dab99e", "#d3ab8b", "#cc9e79", "#be8354", "#b47645")
  } else if(trait == 'Monocyte'){
    c <- "#CC79A7"
    val <- c("#e1b0cb", "#d38bb3", "#CC79A7", "#be548f", "#b44582")
  } else if(trait == 'Lymphocyte'){
    c <- "#009E73"
    val <- c("#00d198", "#00b886", "#009E73", "#006b4e", "#00523b")
  } else if(trait == 'Eosinophil'){
    c <- "#79ccc7"
    val <- c("#9edad6", "#8bd3cf", "#79ccc7", "#54beb8", "#45b4ad")
  } else if(trait == 'MCV'){
    c <- "#79cc9e"
    val <- c("#9edab9", "#8bd3ab", "#79cc9e", "#54be83", "#45b476")
  } else if(trait == 'MCH'){
    c <- "#9e002b"
    val <- c("#eb0040", "#d10039", "#9e002b", "#6b001d", "#520016")
  } else if(trait == 'RBC'){
    c <- "#b20072"
    val <- c("#ff00a3", "#e50093", "#b20072", "#7f0051", "#660041")
  } else if(trait == 'Cystatin_C'){
    c <- "#a7cc79"
    val <- c("#bfda9e", "#b3d38b", "#a7cc79", "#8fbe54", "#82b445")
  } else if(trait == 'Platelet'){
    c <- "#a640b0"
    val <- c("#c371cb", "#bb5ec5", "#a640b0", "#83328b", "#712c78")
  } else if(trait == 'Triglycerides'){
    c <- "#56B4E9"
    val <- c("#9ad2f2", "#83c8ef", "#56B4E9", "#29a0e3", "#1c93d7")
  }
  
  plot <- ggplot(mean_df, aes(x = tile, y = mean_pred)) +
    geom_errorbar(aes(ymin = mean_pred - se, ymax = mean_pred + se, color = tile),
                  width = 0.125, linewidth = 1) +
    geom_point(aes(color = tile), size = 3) +
    geom_line(aes(group = phenotype), size = 0.5, color = c) +
    scale_color_manual(values = val) +
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica"),
          axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
          legend.position = "none") +
    scale_y_reverse() +
    coord_cartesian(xlim = c(1.5125, 4.4875))
  
  return(plot)
}

# Make the plots
plot_pred_income_height <- plot_pred_income(ind_pgs_df, "Height")
plot_pred_income_weight <- plot_pred_income(ind_pgs_df, "Weight")
plot_pred_income_ldl <- plot_pred_income(ind_pgs_df, "LDL")
plot_pred_income_wbc <- plot_pred_income(ind_pgs_df, "WBC")

plot_pred_income_bmi <- plot_pred_income(ind_pgs_df, "BMI")
plot_pred_income_body_fat_perc <- plot_pred_income(ind_pgs_df, "Body_Fat_Perc")
plot_pred_income_monocyte <- plot_pred_income(ind_pgs_df, "Monocyte")
plot_pred_income_lymphocyte <- plot_pred_income(ind_pgs_df, "Lymphocyte")
plot_pred_income_eosinophil <- plot_pred_income(ind_pgs_df, "Eosinophil")
plot_pred_income_mcv <- plot_pred_income(ind_pgs_df, "MCV")
plot_pred_income_mch <- plot_pred_income(ind_pgs_df, "MCH")
plot_pred_income_rbc <- plot_pred_income(ind_pgs_df, "RBC")
plot_pred_income_cystatin_c <- plot_pred_income(ind_pgs_df, "Cystatin_C")
plot_pred_income_platelet <- plot_pred_income(ind_pgs_df, "Platelet")
plot_pred_income_triglycerides <- plot_pred_income(ind_pgs_df, "Triglycerides")

# r2
pheno <- c("Height","Cystatin_C", "Platelet","MCV", 
           "Weight", "MCH", "BMI", "RBC", 
           "Body_Fat_Perc", "Monocyte", "Triglycerides", "Lymphocyte",
           "WBC","Eosinophil", "LDL")

r2 <- read_csv("data/pgs_pred/ind_pred_r2.csv")
r2 <- r2 %>% filter(Method == "Spline" & Predictor %in% c("gen_dist", "townsend"))
r2$Predictor <- factor(r2$Predictor, levels = c("gen_dist", "townsend"), 
                      labels = c("Genetic distance", "Townsend\ndeprivation index"))
r2$Phenotype <- factor(r2$Phenotype, levels = pheno,
                      labels = str_replace_all(pheno, "_", " "))
plot_r2 <- r2 %>% ggplot(aes(x = Predictor, y = R2, group = Phenotype, color = Phenotype)) +
  geom_line(size = 2) +
  geom_point(size = 4) +
  scale_color_manual(values = c("#E69F00", "#a7cc79", "#a640b0", "#79cc9e",
                                "#40B0A6", "#9e002b", "#9e79cc", "#b20072",
                                "#cc9e79", "#CC79A7", "#56B4E9", "#009E73",
                                "#0072B2", "#79ccc7", "#DC3220")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size=24, family = "Helvetica"),
        axis.title.x = element_blank(),
        axis.text = element_text(size=20, family = "Helvetica"),
        axis.text.x = element_text(size=24, family = "Helvetica", color = "black"),
        axis.text.y = element_text(vjust = 0.25, hjust = 0.5),
        legend.position = "none") +
  scale_y_log10(
    breaks = scales::trans_breaks("log10", function(x) 10^x),
    labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(sides = "l", outside = TRUE) +
  coord_cartesian(clip = "off") +
  annotate("text", label = "Height", x = 0.95, y = 10^-2.49, size = 8,  family = "Helvetica",
           color = "#E69F00", hjust = 1, fontface = "bold") +
  annotate("text", label = "Cystatin C", x = 2.05, y = 10^-2.49, size = 8,  family = "Helvetica",
           color = "#a7cc79", hjust = 0, fontface = "bold") +
  annotate("text", label = "Platelet", x = 0.95, y = 10^-2.84, size = 8,  family = "Helvetica",
           color = "#a640b0", hjust = 1, fontface = "bold") +
  annotate("text", label = "MCV", x = 0.95, y = 10^-1.83, size = 8,  family = "Helvetica",
           color = "#79cc9e", hjust = 1, fontface = "bold") +
  annotate("text", label = "Weight", x = 2.05, y = 10^-2.39, size = 8,  family = "Helvetica",
           color = "#40B0A6", hjust = 0, fontface = "bold") +
  annotate("text", label = "MCH", x = 2.05, y = 10^-2.59, size = 8,  family = "Helvetica",
           color = "#9e002b", hjust = 0, fontface = "bold") +
  annotate("text", label = "BMI", x = 2.05, y = 10^-2.31, size = 8,  family = "Helvetica",
           color = "#9e79cc", hjust = 0, fontface = "bold") +
  annotate("text", label = "RBC", x = 0.95, y = 10^-2.23, size = 8,  family = "Helvetica",
           color = "#b20072", hjust = 1, fontface = "bold") +
  annotate("text", label = "Body Fat Perc", x = 2.05, y = 10^-2.23, size = 8,  family = "Helvetica",
           color = "#cc9e79", hjust = 0, fontface = "bold") +
  annotate("text", label = "Monocyte", x = 2.05, y = 10^-2.85, size = 8,  family = "Helvetica",
           color = "#CC79A7", hjust = 0, fontface = "bold") +
  annotate("text", label = "Triglycerides", x = 0.95, y = 10^-2.76, size = 8,  family = "Helvetica",
           color = "#56B4E9", hjust = 1, fontface = "bold") +
  annotate("text", label = "Lymphocyte", x = 2.05, y = 10^-3.0, size = 8,  family = "Helvetica",
           color = "#009E73", hjust = 0, fontface = "bold") +
  annotate("text", label = "WBC", x = 2.05, y = 10^-2.67, size = 8,  family = "Helvetica",
           color = "#0072B2", hjust = 0, fontface = "bold") +
  annotate("text", label = "Eosinophil", x = 0.95, y = 10^-3.02, size = 8,  family = "Helvetica",
           color = "#79ccc7", hjust = 1, fontface = "bold") +
  annotate("text", label = "LDL", x = 2.05, y = 10^-3.1, size = 8,  family = "Helvetica",
           color = "#DC3220", hjust = 0, fontface = "bold") +
  ylab("Explained variance in squared prediction error")

# Combine the panels
fig_height <- grid.arrange(arrangeGrob(plot_pred_dist_height,
                                      top = textGrob("\tHeight",
                                                     gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#E69F00", fontface = "bold"))))

fig_weight <- grid.arrange(arrangeGrob(plot_pred_dist_weight,
                                      top = textGrob("\tWeight",
                                                     gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#40B0A6", fontface = "bold"))))

fig_ldl <- grid.arrange(arrangeGrob(plot_pred_dist_ldl,
                                   top = textGrob("LDL cholesterol level",
                                                  gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#DC3220", fontface = "bold"))))

fig_dist <- plot_grid(fig_height, 
                     fig_weight, 
                     fig_ldl, 
                     #fig_wbc, 
                     ncol = 3, nrow = 1)

fig_townsend <- plot_grid(plot_pred_townsend_height,
                         plot_pred_townsend_weight,
                         plot_pred_townsend_ldl,
                         ncol = 3, nrow = 1)

fig_income <- plot_grid(plot_pred_income_height,
                       plot_pred_income_weight,
                       plot_pred_income_ldl,
                       ncol = 3, nrow = 1)

a <- grid.arrange(arrangeGrob(fig_dist,
                             bottom = textGrob("Genetic distance", 
                                               gp=gpar(fontfamily = "Helvetica", fontsize=24)),
                             top = textGrob("↑ Better prediction\t\t\t\t\t\t\t\t\t\t\t\t",
                                            gp=gpar(fontfamily = "Helvetica", fontsize=24))))
b <- grid.arrange(arrangeGrob(fig_townsend,
                              bottom = textGrob("Townsend deprivation index",
                                                gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                              top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))
c <- grid.arrange(arrangeGrob(fig_income,
                              bottom = textGrob("Household income (£)",
                                                gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                              top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))

fig_3_a <- plot_grid(a, NULL, b, NULL, c,
                    label_x = 0.01, hjust = 0,
                    label_size = 28, scale = 1,
                    rel_heights = c(1.4, 0.1, 1.3, 0.1, 1.6),
                    label_fontfamily = "Helvetica",
                    nrow = 5)

fig_3_a <- grid.arrange(arrangeGrob(fig_3_a,
                                   left = textGrob("Mean squared prediction error", rot = 90,
                                                   gp=gpar(fontfamily = "Helvetica", fontsize=24))))

fig_3 <- plot_grid(NULL, NULL, NULL,
                  fig_3_a, NULL, plot_r2,
                  labels = c('A. Mean trends in individual-level prediction accuracy\n', 
                             '',
                             'B. Deprivation index and genetic distance explain portability\ncomparably well', '', ''), 
                  label_x = 0.01, hjust = 0,
                  label_size = 28, scale = 1,
                  rel_heights = c(0.15, 1),
                  rel_widths = c(1, 0.1, 1),
                  label_fontfamily = "Helvetica",
                  nrow = 2, ncol = 3)

grDevices::cairo_pdf("img/fig_3_ind_pred.pdf", width = 24, height = 12, onefile = T)
print(fig_3)
dev.off()

# Supp fig for mean trends
fig_bmi <- grid.arrange(arrangeGrob(plot_pred_dist_bmi,
                                   top = textGrob("\tBMI",
                                                  gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#9e79cc", fontface = "bold"))))

fig_body_fat_perc <- grid.arrange(arrangeGrob(plot_pred_dist_body_fat_perc,
                                             top = textGrob("Body fat percentage",
                                                            gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#cc9e79", fontface = "bold"))))

fig_monocyte <- grid.arrange(arrangeGrob(plot_pred_dist_monocyte,
                                        top = textGrob("\tMonocyte count",
                                                       gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#CC79A7", fontface = "bold"))))

fig_lymphocyte <- grid.arrange(arrangeGrob(plot_pred_dist_lymphocyte,
                                          top = textGrob("\tLymphocyte count",
                                                         gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#009E73", fontface = "bold"))))

fig_wbc <- grid.arrange(arrangeGrob(plot_pred_dist_wbc,
                                   top = textGrob("White blood cell count",
                                                  gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#0072B2", fontface = "bold"))))

fig_eosinophil <- grid.arrange(arrangeGrob(plot_pred_dist_eosinophil,
                                          top = textGrob("\tEosinphil count",
                                                         gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#79ccc7", fontface = "bold"))))

fig_mcv <- grid.arrange(arrangeGrob(plot_pred_dist_mcv,
                                   top = textGrob("\tMean corpuscular\n\tvolume",
                                                  gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#79cc9e", fontface = "bold"))))

fig_mch <- grid.arrange(arrangeGrob(plot_pred_dist_mch,
                                   top = textGrob("\tMean corpuscular\n\themoglobin",
                                                  gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#9e002b", fontface = "bold"))))

fig_rbc <- grid.arrange(arrangeGrob(plot_pred_dist_rbc,
                                   top = textGrob("\tRed blood cell\n\tcount",
                                                  gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#b20072", fontface = "bold"))))

fig_cystatin_c <- grid.arrange(arrangeGrob(plot_pred_dist_cystatin_c,
                                          top = textGrob("\tCystatin C level",
                                                         gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#a7cc79", fontface = "bold"))))

fig_platelet <- grid.arrange(arrangeGrob(plot_pred_dist_platelet,
                                        top = textGrob("\tPlatelet count",
                                                       gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#a640b0", fontface = "bold"))))

fig_triglycerides <- grid.arrange(arrangeGrob(plot_pred_dist_triglycerides,
                                             top = textGrob("Triglyceride level",
                                                            gp=gpar(fontfamily = "Helvetica", fontsize=24, col = "#56B4E9", fontface = "bold"))))

fig_dist_1 <- plot_grid(fig_bmi, 
                       fig_body_fat_perc, 
                       ncol = 2, nrow = 1)

fig_townsend_1 <- plot_grid(plot_pred_townsend_bmi,
                           plot_pred_townsend_body_fat_perc,
                           ncol = 2, nrow = 1)

fig_income_1 <- plot_grid(plot_pred_income_bmi,
                         plot_pred_income_body_fat_perc,
                         ncol = 2, nrow = 1)

a_1 <- grid.arrange(arrangeGrob(fig_dist_1,
                               bottom = textGrob("Genetic distance", 
                                                 gp=gpar(fontfamily = "Helvetica", fontsize=24))))

b_1 <- grid.arrange(arrangeGrob(fig_townsend_1,
                                bottom = textGrob("Townsend deprivation index",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))
c_1 <- grid.arrange(arrangeGrob(fig_income_1,
                                bottom = textGrob("Household income (£)",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))

fig_supp_mean_1 <- plot_grid(a_1, NULL, b_1, NULL, c_1,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            rel_heights = c(1.3, 0.1, 1.3, 0.1, 1.6),
                            label_fontfamily = "Helvetica",
                            nrow = 5)
fig_supp_mean_1 <- plot_grid(fig_supp_mean_1, NULL,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            label_fontfamily = "Helvetica",
                            ncol = 2)


grDevices::cairo_pdf("img/fig_s14_mean_trend_physical.pdf", width = 16, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_supp_mean_1,
                         left = textGrob("Mean squared prediction error", rot = 90,
                                         gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

fig_dist_2 <- plot_grid(fig_monocyte, 
                       fig_lymphocyte, 
                       fig_wbc, 
                       fig_eosinophil,
                       ncol = 4, nrow = 1)

fig_townsend_2 <- plot_grid(plot_pred_townsend_monocyte,
                           plot_pred_townsend_lymphocyte,
                           plot_pred_townsend_wbc,
                           plot_pred_townsend_eosinophil,
                           ncol = 4, nrow = 1)

fig_income_2 <- plot_grid(plot_pred_income_monocyte,
                         plot_pred_income_lymphocyte,
                         plot_pred_income_wbc,
                         plot_pred_income_eosinophil,
                         ncol = 4, nrow = 1)

a_2 <- grid.arrange(arrangeGrob(fig_dist_2,
                               bottom = textGrob("Genetic distance", 
                                                 gp=gpar(fontfamily = "Helvetica", fontsize=24))))

b_2 <- grid.arrange(arrangeGrob(fig_townsend_2,
                                bottom = textGrob("Townsend deprivation index",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))
c_2 <- grid.arrange(arrangeGrob(fig_income_2,
                                bottom = textGrob("Household income (£)",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))

fig_supp_mean_2 <- plot_grid(a_2, NULL, b_2, NULL, c_2,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            rel_heights = c(1.3, 0.1, 1.3, 0.1, 1.6),
                            label_fontfamily = "Helvetica",
                            nrow = 5)


grDevices::cairo_pdf("img/fig_s15_mean_trend_wbc.pdf", width = 16, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_supp_mean_2,
                         left = textGrob("Mean squared prediction error", rot = 90,
                                         gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

fig_dist_3 <- plot_grid(fig_mcv, 
                       fig_mch, 
                       fig_rbc,
                       ncol = 3, nrow = 1)

fig_townsend_3 <- plot_grid(plot_pred_townsend_mcv,
                           plot_pred_townsend_mch,
                           plot_pred_townsend_rbc,
                           ncol = 3, nrow = 1)

fig_income_3 <- plot_grid(plot_pred_income_mcv,
                         plot_pred_income_mch,
                         plot_pred_income_rbc,
                         ncol = 3, nrow = 1)

a_3 <- grid.arrange(arrangeGrob(fig_dist_3,
                               bottom = textGrob("Genetic distance", 
                                                 gp=gpar(fontfamily = "Helvetica", fontsize=24))))

b_3 <- grid.arrange(arrangeGrob(fig_townsend_3,
                                bottom = textGrob("Townsend deprivation index",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))
c_3 <- grid.arrange(arrangeGrob(fig_income_3,
                                bottom = textGrob("Household income (£)",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))

fig_supp_mean_3 <- plot_grid(a_3, NULL, b_3, NULL, c_3,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            rel_heights = c(1.4, 0.1, 1.3, 0.1, 1.6),
                            label_fontfamily = "Helvetica",
                            nrow = 5)
fig_supp_mean_3 <- plot_grid(fig_supp_mean_3, NULL,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            label_fontfamily = "Helvetica",
                            ncol = 2,
                            rel_widths = c(3, 1))


grDevices::cairo_pdf("img/fig_s16_mean_trend_rbc.pdf", width = 16, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_supp_mean_3,
                         left = textGrob("Mean squared prediction error", rot = 90,
                                         gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

fig_dist_4 <- plot_grid(fig_cystatin_c, 
                       fig_platelet, 
                       fig_triglycerides,
                       ncol = 3, nrow = 1)

fig_townsend_4 <- plot_grid(plot_pred_townsend_cystatin_c,
                           plot_pred_townsend_platelet,
                           plot_pred_townsend_triglycerides,
                           ncol = 3, nrow = 1)

fig_income_4 <- plot_grid(plot_pred_income_cystatin_c,
                         plot_pred_income_platelet,
                         plot_pred_income_triglycerides,
                         ncol = 3, nrow = 1)

a_4 <- grid.arrange(arrangeGrob(fig_dist_4,
                               bottom = textGrob("Genetic distance", 
                                                 gp=gpar(fontfamily = "Helvetica", fontsize=24))))

b_4 <- grid.arrange(arrangeGrob(fig_townsend_4,
                                bottom = textGrob("Townsend deprivation index",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))
c_4 <- grid.arrange(arrangeGrob(fig_income_4,
                                bottom = textGrob("Household income (£)",
                                                  gp = gpar(fontfamily = "Helvetica", fontsize = 24)),
                                top = linesGrob(x = c(0, 1), y = 0.5, gp = gpar(col = "black", lwd = 2))))

fig_supp_mean_4 <- plot_grid(a_4, NULL, b_4, NULL, c_4,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            rel_heights = c(1.3, 0.1, 1.3, 0.1, 1.6),
                            label_fontfamily = "Helvetica",
                            nrow = 5)
fig_supp_mean_4 <- plot_grid(fig_supp_mean_4, NULL,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            label_fontfamily = "Helvetica",
                            ncol = 2,
                            rel_widths = c(3, 1))


grDevices::cairo_pdf("img/fig_s17_mean_trend_other.pdf", width = 16, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_supp_mean_4,
                         left = textGrob("Mean squared prediction error", rot = 90,
                                         gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


# Supp fig for r2
r2 <- read_csv("data/pgs_pred/ind_pred_r2.csv")
r2 <- r2 %>% separate(`95% CI`, c("lower", "upper"), ", ")
r2 <- r2 %>% mutate(lower = as.numeric(str_remove(lower, "\\(")),
                   upper = as.numeric(str_remove(upper, "\\)")))
r2$Predictor <- factor(r2$Predictor, levels = c("townsend", "year_edu", "income", "gen_dist", 
                                               "ma_count_low", "ma_count_medium", "ma_count_high",
                                               "ma_count_total"), 
                      labels = c("Townsend\ndeprivation\nindex", 
                                 "Years of\neducation", "Household\nincome (£)",
                                 "Genetic distance", 
                                 "Minor allele count\nof small effect\nsize SNPs",
                                 "Minor allele count\nof medium effect\nsize SNPs",
                                 "Minor allele count\nof large effect\nsize SNPs",
                                 "Minor allele count\nof all SNPs"))
r2$Method <- factor(r2$Method, labels = c("Discretized (uniform)", "Linear", "Spline"))
# Plot as %
r2$R2 <- r2$R2 * 100
r2$lower <- r2$lower * 100
r2$upper <- r2$upper * 100

plot_r2_supp <- function(r2, trait){
  r2 <- r2 %>% filter(Phenotype == trait)
  
  plot <- r2 %>% ggplot(aes(x = Method, y = R2, color = Predictor)) +
    geom_hline(yintercept = 0, size = 0.5, color = "black") +
    geom_point(size = 2.5, position = position_dodge(width = 1)) +
    geom_errorbar(aes(ymin = lower, ymax = upper),
                  width = 0.5, linewidth = 1, position = position_dodge(width = 1)) +
    geom_vline(xintercept = 1.5, linetype = "dashed", color = "black", linewidth = 0.5) +
    geom_vline(xintercept = 2.5, linetype = "dashed", color = "black", linewidth = 0.5) +
    theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(size=24, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica"),
          axis.text.x = element_text(size=24, family = "Helvetica", color = "black"),
          axis.text.y = element_text(vjust = 0.25, hjust = 0.5),
          legend.position = "bottom",
          legend.text = element_text(size=18, family = "Helvetica"),
          legend.title = element_text(size=20, family = "Helvetica")) +
    scale_color_manual(values = c("#4b132f", "#AA4499", "#CC6677", "#DDCC77", "#88CCEE",
                                  "#44AA99", '#117733', "#332288")) +
    scale_y_continuous(labels = label_number()) +
    ylab("Fraction of variance\nexplained in squared\nprediction error (%)")
  
  return(plot)
}

# Make the plots
plot_r2_supp_height <- plot_r2_supp(r2, "Height")
plot_r2_supp_weight <- plot_r2_supp(r2, "Weight")
plot_r2_supp_bmi <- plot_r2_supp(r2, "BMI")
plot_r2_supp_body_fat_perc <- plot_r2_supp(r2, "Body_Fat_Perc")

plot_r2_supp_monocyte <- plot_r2_supp(r2, "Monocyte")
plot_r2_supp_lymphocyte <- plot_r2_supp(r2, "Lymphocyte")
plot_r2_supp_wbc <- plot_r2_supp(r2, "WBC")
plot_r2_supp_eosinophil <- plot_r2_supp(r2, "Eosinophil")

plot_r2_supp_mcv <- plot_r2_supp(r2, "MCV")
plot_r2_supp_mch <- plot_r2_supp(r2, "MCH")
plot_r2_supp_rbc <- plot_r2_supp(r2, "RBC")

plot_r2_supp_cystatin_c <- plot_r2_supp(r2, "Cystatin_C")
plot_r2_supp_platelet <- plot_r2_supp(r2, "Platelet")
plot_r2_supp_triglycerides <- plot_r2_supp(r2, "Triglycerides")
plot_r2_supp_ldl <- plot_r2_supp(r2, "LDL")

fig_r2_supp_1 <- plot_grid(NULL, NULL, plot_r2_supp_height, plot_r2_supp_weight, 
                           NULL, NULL, plot_r2_supp_bmi, plot_r2_supp_body_fat_perc, 
                           labels = c('A. Height', 
                                      'B. Weight', 
                                      '',
                                      '',
                                      'C. BMI', 
                                      'D. Body fat percentage', 
                                      '',
                                      ''), ncol = 2, nrow = 4,
                           label_x = 0.01, hjust = 0,
                           label_size = 28, scale = 1,
                           rel_heights = c(0.15, 1, 0.15, 1),
                           label_fontfamily = "Helvetica")

fig_r2_supp_2 <- plot_grid(NULL, NULL, plot_r2_supp_monocyte, plot_r2_supp_lymphocyte, 
                           NULL, NULL, plot_r2_supp_wbc, plot_r2_supp_eosinophil, 
                           labels = c('A. Monocyte count', 
                                      'B. Lymphocyte count', 
                                      '',
                                      '',
                                      'C. White blood cell count', 
                                      'D. Eosinophil count', 
                                      '',
                                      ''), ncol = 2, nrow = 4,
                           label_x = 0.01, hjust = 0,
                           label_size = 28, scale = 1,
                           rel_heights = c(0.15, 1, 0.15, 1),
                           label_fontfamily = "Helvetica")

fig_r2_supp_3 <- plot_grid(NULL, NULL, plot_r2_supp_mcv, plot_r2_supp_mch, 
                           NULL, NULL, plot_r2_supp_rbc, NULL, 
                           labels = c('A. Mean corpuscular volume', 
                                      'B. Mean corpuscular hemoglobin', 
                                      '',
                                      '',
                                      'C. Red blood cell count', 
                                      '', 
                                      '',
                                      ''), ncol = 2, nrow = 4,
                           label_x = 0.01, hjust = 0,
                           label_size = 28, scale = 1,
                           rel_heights = c(0.15, 1, 0.15, 1),
                           label_fontfamily = "Helvetica")

fig_r2_supp_4 <- plot_grid(NULL, NULL, plot_r2_supp_cystatin_c, plot_r2_supp_platelet, 
                           NULL, NULL, plot_r2_supp_triglycerides, plot_r2_supp_ldl, 
                           labels = c('A. Cystatin C level', 
                                      'B. Platelet count', 
                                      '',
                                      '',
                                      'C. Triglyceride level', 
                                      'D. LDL cholesterol level', 
                                      '',
                                      ''), ncol = 2, nrow = 4,
                           label_x = 0.01, hjust = 0,
                           label_size = 28, scale = 1,
                           rel_heights = c(0.15, 1, 0.15, 1),
                           label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s18_predictor_r2_physical.pdf", width = 24, height = 12, onefile = T)
print(fig_r2_supp_1)
dev.off()

grDevices::cairo_pdf("img/fig_s19_predictor_r2_wbc.pdf", width = 24, height = 12, onefile = T)
print(fig_r2_supp_2)
dev.off()

grDevices::cairo_pdf("img/fig_s20_predictor_r2_rbc.pdf", width = 24, height = 12, onefile = T)
print(fig_r2_supp_3)
dev.off()

grDevices::cairo_pdf("img/fig_s21_predictor_r2_other.pdf", width = 24, height = 12, onefile = T)
print(fig_r2_supp_4)
dev.off()
