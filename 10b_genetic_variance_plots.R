library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)

`%notin%` <- Negate(`%in%`)

pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
          "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
          "Body_Fat_Perc")

# Effect sizes in GWAS set
snp_lst <- list()

for(trait in pheno){
  snp <- read_tsv(paste0("/work/06568/joyce_w/stampede2/pgs_portability/data/prs/", trait, 
                        "_imputation_array_WB_filter_no_std_2nd_round_threshold_1.txt"), 
                 col_names = F)
  colnames(snp) <- "ID"
  snp <- snp[snp$ID %notin% mhc,]
  
  beta <- read_tsv(paste0("/work/06568/joyce_w/stampede2/pgs_portability/data/gwas_results/", trait, 
                         "_imputation_array_WB_filter_no_std_2nd_round_combined.glm.linear"))
  snp <- snp %>% left_join(beta, by = "ID")
  snp$beta_sq <- (snp$BETA)^2
  
  snp <- snp[, c(1, 14)]
  
  snp_lst[[trait]] <- snp
}

snp_beta <- snp_lst[[pheno[1]]]
snp_beta$category <- ntile(snp_beta$beta_sq, 3)
snp_beta$phenotype <- pheno[1]

for(trait in pheno[2:15]){
  snp_beta_temp <- snp_lst[[trait]]
  snp_beta_temp$category <- ntile(snp_beta_temp$beta_sq, 3)
  snp_beta_temp$phenotype <- trait
  snp_beta <- snp_beta %>% rbind.data.frame(snp_beta_temp)
}

# Phynotypic variance

gwas_id <- read_table('data/ukb_populations/wb_gwas_id.txt')

phenotype <- read_table("data/phenotypes/phenotypes.tsv")
phenotype <- phenotype %>% select(-1)
phenotype <- phenotype %>% filter(IID %in% gwas_id$IID) %>%
  pivot_longer(cols = BMI:Body_Fat_Perc, names_to = "phenotype", values_to = "phenotype_value")

phenotype_var <- phenotype %>% group_by(phenotype) %>%
  summarize(phenotype_var = var(phenotype_value, na.rm = T))

snp_beta <- snp_beta %>% left_join(phenotype_var, by = "phenotype")
snp_beta$beta_sq_std <- snp_beta$beta_sq / snp_beta$phenotype_var

# Read close and far GWAS betas
close <- data.frame()
far <- data.frame()
for(trait in pheno){
  close_temp <- read_tsv(paste0("data/gwas_results/", trait, "_close.", trait, ".glm.linear"))
  far_temp <- read_tsv(paste0("data/gwas_results/", trait, "_far.", trait, ".glm.linear"))
  close_temp$beta_sq_close <- (close_temp$BETA) ^ 2
  close_temp <- close_temp %>% select(ID, beta_sq_close)
  close_temp$phenotype <- trait
  
  far_temp$beta_sq_far <- (far_temp$BETA) ^ 2
  far_temp <- far_temp %>% select(ID, beta_sq_far)
  far_temp$phenotype <- trait
  
  close <- rbind.data.frame(close, close_temp)
  far <- rbind.data.frame(far, far_temp)
}

snp_beta <- snp_beta %>% left_join(close, by = c("phenotype", "ID"))
snp_beta <- snp_beta %>% left_join(far, by = c("phenotype", "ID"))


# Calculate heterozygosity
# Heterozygosity (expected) in GWAS
heterozygosity <- read_tsv("data/heterozygosity/heterozygosity_wb.frqx")
heterozygosity$total <- heterozygosity$`C(HOM A1)` + heterozygosity$`C(HET)` + heterozygosity$`C(HOM A2)`
heterozygosity$p <- heterozygosity$`C(HOM A1)` / heterozygosity$total + (heterozygosity$`C(HET)`) / 2 / heterozygosity$total
heterozygosity$heterozygosity <- 2 * heterozygosity$p * (1 - heterozygosity$p)
heterozygosity <- heterozygosity %>% select(SNP, heterozygosity)
colnames(heterozygosity)[1] <- "ID"
heterozygosity %>% write_tsv("data/heterozygosity/heterozygosity_wb.tsv")

colnames(heterozygosity)[2] <- "het_gwas"
snp_beta <- snp_beta %>% left_join(heterozygosity, by = "ID")
snp_beta$beta_sq_by_het_gwas <- snp_beta$beta_sq_std * snp_beta$het_gwas

# Heterozygosity (expected) in close
heterozygosity <- read_tsv("data/heterozygosity/heterozygosity_close.frqx")
heterozygosity$total <- heterozygosity$`C(HOM A1)` + heterozygosity$`C(HET)` + heterozygosity$`C(HOM A2)`
heterozygosity$p <- heterozygosity$`C(HOM A1)` / heterozygosity$total + (heterozygosity$`C(HET)`) / 2 / heterozygosity$total
heterozygosity$heterozygosity <- 2 * heterozygosity$p * (1 - heterozygosity$p)
heterozygosity <- heterozygosity %>% select(SNP, heterozygosity)
colnames(heterozygosity)[1] <- "ID"
heterozygosity %>% write_tsv("data/heterozygosity/heterozygosity_close.tsv")

colnames(heterozygosity)[2] <- "het_close"
snp_beta <- snp_beta %>% left_join(heterozygosity, by = "ID")
snp_beta$beta_sq_by_het_close <- snp_beta$beta_sq_close / snp_beta$phenotype_var * snp_beta$het_close

# Heterozygosity (expected) in far
heterozygosity <- read_tsv("data/heterozygosity/heterozygosity_far.frqx")
heterozygosity$total <- heterozygosity$`C(HOM A1)` + heterozygosity$`C(HET)` + heterozygosity$`C(HOM A2)`
heterozygosity$p <- heterozygosity$`C(HOM A1)` / heterozygosity$total + (heterozygosity$`C(HET)`) / 2 / heterozygosity$total
heterozygosity$heterozygosity <- 2 * heterozygosity$p * (1 - heterozygosity$p)
heterozygosity <- heterozygosity %>% select(SNP, heterozygosity)
colnames(heterozygosity)[1] <- "ID"
heterozygosity %>% write_tsv("data/heterozygosity/heterozygosity_far.tsv")

colnames(heterozygosity)[2] <- "het_far"
snp_beta <- snp_beta %>% left_join(heterozygosity, by = "ID")
snp_beta$beta_sq_by_het_far <- snp_beta$beta_sq_far / snp_beta$phenotype_var * snp_beta$het_far

snp_beta$category <- factor(snp_beta$category, labels = c("small", "medium", "large"))
snp_beta$phenotype <- factor(snp_beta$phenotype, 
                            levels = rev(c("Height", "Platelet", "MCV", "MCH", "RBC", 
                                           "Monocyte", "LDL", "Triglycerides", "Cystatin_C", "Weight",
                                           "WBC", "Eosinophil", "BMI", "Body_Fat_Perc", "Lymphocyte")),
                            labels = rev(c("Height (h\u00b2 = 0.49)", 
                                           "Platelet (h\u00b2 = 0.31)", 
                                           "MCV (h\u00b2 = 0.27)", 
                                           "MCH (h\u00b2 = 0.25)", 
                                           "RBC (h\u00b2 = 0.23)", 
                                           "Monocyte (h\u00b2 = 0.23)",
                                           "LDL (h\u00b2 = 0.08)",
                                           "Triglycerides (h\u00b2 = 0.22)", 
                                           "Cystatin C (h\u00b2 = 0.32)", 
                                           "Weight (h\u00b2 = 0.27)", 
                                           "WBC (h\u00b2 = 0.19)", 
                                           "Eosinophil (h\u00b2 = 0.18)",
                                           "BMI (h\u00b2 = 0.25)", 
                                           "Body fat perc (h\u00b2 = 0.23)", 
                                           "Lymphocyte (h\u00b2 = 0.21)")))

# Make plots
plot_beta_sq <- snp_beta %>% 
  ggplot(aes(x = beta_sq_std, y = phenotype, color = category)) + 
  geom_point(size = 2, alpha = 0.05, position = position_dodge(width=0.9)) +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        axis.text.y=element_text(size=20, family = "Helvetica", color = "black"),
        plot.title=element_text(size=28, face = "bold", family = "Helvetica"),
        plot.subtitle=element_text(size=18, face = "bold", family = "Helvetica"),
        plot.caption=element_text(size=16, family = "Helvetica"),
        axis.title.y=element_blank(),
        legend.position = "none") +
  xlab("Squared allelic effect estimate") +
  guides(color = guide_legend(title="Effect size sq.")) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(base = 10, sides = "b", outside = T) + 
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c("#ffcd34", "#ff8934", "#ff3434")) +
  annotate("text", label = "Effect size", x = 10^-0.1, y = 3.4, size = 8,  family = "Helvetica",
           color = "black", hjust = 0) +
  annotate("text", label = "small", x = 10^-0.1, y = 2.6, size = 8,  family = "Helvetica",
           color = "#ffcd34", hjust = 0) +
  annotate("text", label = "medium", x = 10^-0.1, y = 1.8, size = 8,  family = "Helvetica",
           color = "#ff8934", hjust = 0) +
  annotate("text", label = "large", x = 10^-0.1, y = 1, size = 8,  family = "Helvetica",
           color = "#ff3434", hjust = 0)

pdf(file = "img/fig_s18_beta_sq.pdf", width = 12, height = 6)
print(plot_beta_sq)
dev.off()

snp_beta_2 <- snp_beta %>%
  pivot_longer(cols = c(beta_sq_by_het_close, beta_sq_by_het_far, beta_sq_by_het_gwas), 
               names_to = "group", values_to = "beta_sq_by_het")
snp_beta_2$group <- ifelse(snp_beta_2$group == "beta_sq_by_het_close", "close", 
                          ifelse(snp_beta_2$group == "beta_sq_by_het_far", "far","GWAS"))
snp_beta_2$group <- factor(snp_beta_2$group, levels = c("GWAS", "close", "far"))
snp_beta_2$labels <- paste(snp_beta_2$phenotype, snp_beta_2$group, sep = " in ")
snp_beta_2$labels <- factor(snp_beta_2$labels,
                           levels = rev(c("Height (h² = 0.49) in GWAS", 
                                          "Height (h² = 0.49) in close",
                                          "Height (h² = 0.49) in far",
                                          "Platelet (h² = 0.31) in GWAS",
                                          "Platelet (h² = 0.31) in close",
                                          "Platelet (h² = 0.31) in far",
                                          "MCV (h² = 0.27) in GWAS",
                                          "MCV (h² = 0.27) in close",
                                          "MCV (h² = 0.27) in far",
                                          "MCH (h² = 0.25) in GWAS",
                                          "MCH (h² = 0.25) in close",
                                          "MCH (h² = 0.25) in far",
                                          "RBC (h² = 0.23) in GWAS",
                                          "RBC (h² = 0.23) in close",
                                          "RBC (h² = 0.23) in far",
                                          "Monocyte (h² = 0.23) in GWAS",
                                          "Monocyte (h² = 0.23) in close",
                                          "Monocyte (h² = 0.23) in far",
                                          "LDL (h² = 0.08) in GWAS",
                                          "LDL (h² = 0.08) in close",
                                          "LDL (h² = 0.08) in far",
                                          "Triglycerides (h² = 0.22) in GWAS",
                                          "Triglycerides (h² = 0.22) in close",
                                          "Triglycerides (h² = 0.22) in far",
                                          "Cystatin C (h² = 0.32) in GWAS",
                                          "Cystatin C (h² = 0.32) in close",
                                          "Cystatin C (h² = 0.32) in far",
                                          "Weight (h² = 0.27) in GWAS",
                                          "Weight (h² = 0.27) in close",
                                          "Weight (h² = 0.27) in far",
                                          "WBC (h² = 0.19) in GWAS",
                                          "WBC (h² = 0.19) in close",
                                          "WBC (h² = 0.19) in far",
                                          "Eosinophil (h² = 0.18) in GWAS",
                                          "Eosinophil (h² = 0.18) in close",
                                          "Eosinophil (h² = 0.18) in far",
                                          "BMI (h² = 0.25) in GWAS",
                                          "BMI (h² = 0.25) in close",
                                          "BMI (h² = 0.25) in far",
                                          "Body fat perc (h² = 0.23) in GWAS",
                                          "Body fat perc (h² = 0.23) in close",
                                          "Body fat perc (h² = 0.23) in far",
                                          "Lymphocyte (h² = 0.21) in GWAS",
                                          "Lymphocyte (h² = 0.21) in close",
                                          "Lymphocyte (h² = 0.21) in far")),
                           labels = rev(c("Height (h² = 0.49) in GWAS", 
                                          "Height in \"close\"",
                                          "Height in \"far\"",
                                          "Platelet (h² = 0.31) in GWAS",
                                          "Platelet in \"close\"",
                                          "Platelet in \"far\"",
                                          "MCV (h² = 0.27) in GWAS",
                                          "MCV in \"close\"",
                                          "MCV in \"far\"",
                                          "MCH (h² = 0.25) in GWAS",
                                          "MCH in \"close\"",
                                          "MCH in \"far\"",
                                          "RBC (h² = 0.23) in GWAS",
                                          "RBC in \"close\"",
                                          "RBC in \"far\"",
                                          "Monocyte (h² = 0.23) in GWAS",
                                          "Monocyte in \"close\"",
                                          "Monocyte in \"far\"",
                                          "LDL (h² = 0.08) in GWAS",
                                          "LDL in \"close\"",
                                          "LDL in \"far\"",
                                          "Triglycerides (h² = 0.22) in GWAS",
                                          "Triglycerides in \"close\"",
                                          "Triglycerides in \"far\"",
                                          "Cystatin C (h² = 0.32) in GWAS",
                                          "Cystatin C in \"close\"",
                                          "Cystatin C in \"far\"",
                                          "Weight (h² = 0.27) in GWAS",
                                          "Weight in \"close\"",
                                          "Weight in \"far\"",
                                          "WBC (h² = 0.19) in GWAS",
                                          "WBC \"close\"",
                                          "WBC \"far\"",
                                          "Eosinophil (h² = 0.18) in GWAS",
                                          "Eosinophil \"close\"",
                                          "Eosinophil \"far\"",
                                          "BMI (h² = 0.25) in GWAS",
                                          "BMI in \"close\"",
                                          "BMI in \"far\"",
                                          "Body fat perc (h² = 0.23) in GWAS",
                                          "Body fat perc in \"close\"",
                                          "Body fat perc in \"far\"",
                                          "Lymphocyte (h² = 0.21) in GWAS",
                                          "Lymphocyte in \"close\"",
                                          "Lymphocyte in \"far\"")))

plot_genetic_variance <- snp_beta_2 %>% 
  ggplot(aes(x = beta_sq_by_het, y = labels, color = category)) + 
  geom_point(size = 2, alpha = 0.05, position = position_dodge(width=0.9)) +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        axis.text.y=element_text(size=20, family = "Helvetica", color = "black"),
        plot.title=element_text(size=28, face = "bold", family = "Helvetica"),
        plot.subtitle=element_text(size=18, face = "bold", family = "Helvetica"),
        plot.caption=element_text(size=16, family = "Helvetica"),
        axis.title.y=element_blank(),
        legend.position = "none") +
  xlab("Genetic variance") +
  guides(color = guide_legend(title="Effect size sq.")) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(base = 10, sides = "b", outside = T) + 
  coord_cartesian(clip = "off") +
  scale_color_manual(values = c("#ffcd34", "#ff8934", "#ff3434")) +
  annotate("text", label = "Effect size", x = 10^-12, y = 3.4, size = 8,  family = "Helvetica",
           color = "black", hjust = 0) +
  annotate("text", label = "small", x = 10^-12, y = 2.6, size = 8,  family = "Helvetica",
           color = "#ffcd34", hjust = 0) +
  annotate("text", label = "medium", x = 10^-12, y = 1.8, size = 8,  family = "Helvetica",
           color = "#ff8934", hjust = 0) +
  annotate("text", label = "large", x = 10^-12, y = 1, size = 8,  family = "Helvetica",
           color = "#ff3434", hjust = 0)

pdf(file = "img/fig_s19_genetic_variance.pdf", width = 12, height = 18)
print(plot_genetic_variance)
dev.off()