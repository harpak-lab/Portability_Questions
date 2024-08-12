library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(cowplot)
library(grid)
library(gridExtra)
library(ggokabeito)
library(stringr)

pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
          "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
          "Body_Fat_Perc")

# First compare effect sizes
gwas_beta <- list()
mean <- c()
sd <- c()

for(trait in pheno){
  snp <- read.table(paste0("data/pgs/", trait, "_threshold_1.txt"))
  gwas_snp <- read_tsv(paste0("data/gwas_results/", trait, "_combined.glm.linear"))
  gwas_snp <- gwas_snp %>% filter(ID %in% snp$V1)
  gwas_snp <- gwas_snp %>% arrange(ID)
  
  # Make sure the effect alleles are the same
  # If not, then flip the sign
  close_snp <- read_tsv(paste0("data/gwas_results/", trait, "_close.", trait, ".glm.linear"))
  close_snp <- close_snp %>% filter(ID %in% snp$V1)
  close_snp <- close_snp %>% arrange(ID)
  close_snp$effect_allele <- gwas_snp$A1
  close_snp$BETA <- ifelse(close_snp$A1 == close_snp$effect_allele, close_snp$BETA, -close_snp$BETA)
  
  far_snp <- read_tsv(paste0("data/gwas_results/", trait, "_far.", trait, ".glm.linear"))
  far_snp <- far_snp %>% filter(ID %in% snp$V1)
  far_snp <- far_snp %>% arrange(ID)
  far_snp$effect_allele <- gwas_snp$A1
  far_snp$BETA <- ifelse(far_snp$A1 == far_snp$effect_allele, far_snp$BETA, -far_snp$BETA)
  
  beta <- cbind.data.frame(ID = gwas_snp$ID, GWAS = gwas_snp$BETA)
  
  close_snp <- close_snp %>% select(ID, BETA)
  colnames(close_snp) <- c("ID", "close")
  
  far_snp = far_snp %>% select(ID, BETA)
  colnames(far_snp) <- c("ID", "far")
  
  beta <- beta %>% left_join(close_snp, by = "ID")
  beta <- beta %>% left_join(far_snp, by = "ID")
  
  # Calculate the betas relative to the original GWAS
  beta <- na.omit(beta)
  beta$GWAS_rel <- beta$GWAS / beta$GWAS
  beta$close_rel <- beta$close / beta$GWAS
  beta$far_rel <- beta$far / beta$GWAS
  
  # Get mean and SD
  mean <- c(mean, mean(beta$GWAS_rel), mean(beta$close_rel), mean(beta$far_rel))
  sd <- c(sd, sd(beta$GWAS_rel), sd(beta$close_rel), sd(beta$far_rel))
  
  beta <- beta %>% pivot_longer(cols = GWAS_rel:far_rel, names_to = "group", values_to = "beta_rel")
  beta$group <- ifelse(beta$group == "GWAS_rel", "original GWAS",
                      ifelse(beta$group == "close_rel", "close", "far"))
  beta$group <- factor(beta$group, levels = c("original GWAS", "close", "far"))
  
  gwas_beta[[trait]] <- beta
}

mean_sd <- cbind.data.frame(phenotype = rep(pheno, each = 3),
                           group = rep(c("original GWAS", "close", "far"), 15),
                           mean = mean,
                           sd = sd)
mean_sd$group <- factor(mean_sd$group, levels = c("original GWAS", "close", "far"))


beta_plot_df <- rbind.data.frame(cbind.data.frame(gwas_beta$Lymphocyte, phenotype = "Lymphocyte"), 
                                cbind.data.frame(gwas_beta$Height, phenotype = "Height"), 
                                cbind.data.frame(gwas_beta$Triglycerides, phenotype = "Triglycerides"))
beta_plot_df$phenotype <- factor(beta_plot_df$phenotype, 
                                levels = c("Height", "Triglycerides", "Lymphocyte"))

# Plot for comparing effect sizes
plot_compare_beta <- beta_plot_df %>% 
  ggplot(aes(x = group, color = phenotype)) + 
  geom_hline(yintercept = 1) +
  geom_errorbar(data = subset(mean_sd, phenotype %in% c("Lymphocyte", "Height", "Triglycerides")) %>% 
                  mutate(phenotype = factor(phenotype, levels = c("Height", "Triglycerides", "Lymphocyte"))), 
                aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.2, linewidth = 1.5,
                position = position_dodge(width = 0.5)) +
  geom_point(data = subset(mean_sd, phenotype %in% c("Lymphocyte", "Height", "Triglycerides")) %>% 
               mutate(phenotype = factor(phenotype, levels = c("Height", "Triglycerides", "Lymphocyte"))), 
             aes(y = mean), size = 4.5,
             position = position_dodge(width = 0.5)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size=24, family = "Helvetica"),
        axis.text = element_text(size=20, family = "Helvetica"),
        axis.title.x = element_blank(),
        legend.position ="none") +
  ylab("Allelic effect estimate,\nrelative to original GWAS sample") +
  labs(color = "Phenotype") +
  scale_color_okabe_ito() +
  annotate("text", label = "Height", x = 0.7, y = 0.8, size = 8,  family = "Helvetica",
           color = "#E69F00") +
  annotate(geom = "segment", x = 0.75, y = 0.88, xend = 0.8, yend = 0.95, 
           color = "#E69F00", size = 1) +
  annotate("text", label = "Triglycerides", x = 1, y = 1.2, size = 8,  family = "Helvetica",
           color = "#56B4E9") +
  annotate(geom = "segment", x = 1, y = 1.12, xend = 1, yend = 1.05, 
           color = "#56B4E9", size = 1) +
  annotate("text", label = "Lymphocyte", x = 1.4, y = 0.8, size = 8,  family = "Helvetica",
           color = "#009E73") +
  annotate(geom = "segment", x = 1.25, y = 0.88, xend = 1.2, yend = 0.95, 
           color = "#009E73", size = 1)

# Then compare heterozygosity
heterozygosity <- read_tsv("data/heterozygosity/mean_heterozygosity.tsv")
heterozygosity <- heterozygosity %>% 
  filter(phenotype %in% c("Lymphocyte", "Height", "Triglycerides"))
heterozygosity$size <- ifelse(heterozygosity$size == 1, "small",
                             ifelse(heterozygosity$size == 2, "medium", "large"))
heterozygosity$size <- factor(heterozygosity$size, levels = c("small", "medium", "large"))

# Upper limit for the GWAS set
pc_dist_gwas <- read_tsv("data/pca/pc_dist_best_gwas_std.tsv")
upper <- unname(quantile(pc_dist_gwas$pc_dist, 0.975))

# Prepare data frames for plot
heterozygosity <- heterozygosity %>% filter(median_pc > upper)
heterozygosity$phenotype <- factor(heterozygosity$phenotype,
                                  levels = c("Height", "Triglycerides", "Lymphocyte"))
heterozygosity$forcolor <- paste(heterozygosity$phenotype, heterozygosity$size)
heterozygosity$forcolor <- factor(heterozygosity$forcolor, 
                                 levels = c("Height small", "Height medium", "Height large",
                                            "Triglycerides small", "Triglycerides medium", "Triglycerides large",
                                            "Lymphocyte small", "Lymphocyte medium", "Lymphocyte large"))

heterozygosity_label <- cbind.data.frame(phenotype = c(rep("Height", 3), 
                                                      rep("Triglycerides", 3),
                                                      rep("Lymphocyte", 3)),
                                        label = rep(c("small effect",
                                                      "medium effect",
                                                      "large effect"), 3),
                                        forcolor = unique(heterozygosity$forcolor),
                                        median_pc = 200,
                                        mean_het = c(0.41, 0.17, 0.03, 
                                                     0.42, 0.16, 0.03,
                                                     0.39, 0.17, 0))
heterozygosity_label$phenotype <- factor(heterozygosity_label$phenotype,
                                        levels = c("Height", "Triglycerides", "Lymphocyte"))
trait_label <- heterozygosity_label
trait_label$label <- c(rep("Height", 3), rep("Triglycerides", 3), rep("Lymphocyte", 3))
trait_label$median_pc <- 5
trait_label$mean_het <- 0.45

heterozygosity_segment <- heterozygosity_label
heterozygosity_segment$median_pc <- 150
heterozygosity_segment$y <- c(0.39, 0.19, 0.05,
                              0.4, 0.18, 0.05,
                              0.37, 0.15, 0.01)
heterozygosity_segment$yend <- c(0.36, 0.22, 0.08,
                                 0.37, 0.21, 0.08,
                                 0.34, 0.12, 0.02)

# Compare heterozygosity
plot_heterozygosity <- heterozygosity %>%
  ggplot(aes(x = median_pc, y = mean_het, color = forcolor, fill = forcolor)) +
  scale_x_continuous(breaks=c(0, 40, 80, 120, 160),
                     expand = c(0, 0)) +
  geom_point(size=5, alpha=0.4, shape = 23) +
  guides(color = guide_legend(title="Effect size sq.", nrow = 1),
         fill = guide_legend(title="Effect size sq.", nrow = 1)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size=24, family = "Helvetica"),
        axis.text = element_text(size=20, family = "Helvetica"),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.position="none",
        axis.title.x = element_blank()) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Heterozygosity") +
  geom_text(data = heterozygosity_label, aes(label = label), hjust = 1, size = 8,
            color = c("#ffc034", "#E69F00", "#9a6a00", 
                      "#9ad2f2", "#56B4E9", "#1c93d7", 
                      "#00d198", "#009E73", "#00523b"), 
            family = "Helvetica") +
  geom_text(data = trait_label, aes(label = label), 
            hjust = 0, size = 8, color = c("#E69F00", "#E69F00", "#E69F00",
                                           "#56B4E9", "#56B4E9", "#56B4E9",
                                           "#009E73", "#009E73", "#009E73"),
            family = "Helvetica") +
  geom_segment(data = heterozygosity_segment, aes(xend = median_pc, y = y, yend = yend), size = 1,
               color = c("#ffc034", "#E69F00", "#9a6a00", 
                         "#9ad2f2", "#56B4E9", "#1c93d7", 
                         "#00d198", "#009E73", "#00523b")) +
  facet_wrap(~phenotype) +
  scale_color_manual(values = c("#ffc034", "#E69F00", "#9a6a00", 
                                "#9ad2f2", "#56B4E9", "#1c93d7", 
                                "#00d198", "#009E73", "#00523b")) +
  scale_fill_manual(values = c("#ffc034", "#E69F00", "#9a6a00", 
                               "#9ad2f2", "#56B4E9", "#1c93d7", 
                               "#00d198", "#009E73", "#00523b")) +
  coord_cartesian(xlim = c(upper, 200))

# Finally, look at variance of PGS
pgs_df <- read_tsv('data/pgs_pred/ind_pgs_df.tsv')

# Calculate variance of PGS, relative to the 50 reference bins
var_pgs <- pgs_df %>%
  dplyr::group_by(phenotype, weighted_pc_groups, group_close_to_gwas) %>%
  dplyr::summarize(median_pc = median(pc_dist),
                   var_pgs = var(pgs, na.rm = T))

var_pgs_baseline <- var_pgs %>%
  dplyr::filter(group_close_to_gwas <= 50) %>%
  dplyr::group_by(phenotype) %>%
  dplyr::summarize(var_pgs_baseline = mean(var_pgs))

var_pgs <- var_pgs %>%
  left_join(var_pgs_baseline, by = "phenotype") %>%
  dplyr::mutate(var_pgs_rel = var_pgs / var_pgs_baseline) %>%
  select(-var_pgs_baseline)

# Compare variance of PGS
plot_var_pgs <- var_pgs %>%
  mutate(phenotype = factor(phenotype, 
                            levels = c("Height", "Cystatin_C", "Platelet", "MCV",
                                       "Weight", "MCH", "BMI", "RBC", "Body_Fat_Perc",
                                       "Monocyte", "Triglycerides", "Lymphocyte",
                                       "WBC", "Eosinophil", "LDL"))) %>%
  ggplot(aes(x = median_pc, y = var_pgs_rel, color = phenotype, fill = phenotype)) +
  scale_x_continuous(breaks = c(0, 40, 80, 120, 160),
                     expand = c(0, 0)) +
  geom_hline(yintercept = 1, size = 2) +
  # Only add data points for lymphocyte
  geom_point(data = subset(var_pgs, phenotype == "Lymphocyte"), 
             size = 5, alpha = 0.4, shape = 23) +
  geom_line(aes(linewidth = phenotype), method = "lm", se=F, na.rm=T, stat = "smooth") +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title=element_text(size = 24, family = "Helvetica"),
        axis.text=element_text(size = 20, family = "Helvetica"),
        legend.position="none",
        axis.title.x = element_blank()) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Variance in PGS") +
  scale_linewidth_manual(values = c(0.25, 0.25, 0.25, 0.25,
                                    0.25, 0.25, 0.25, 0.25, 0.25,
                                    2.5, 0.25, 2.5,
                                    2.5, 0.25, 0.25)) +
  scale_color_manual(values = c("#BEBEBE", "#BEBEBE", "#BEBEBE", "#BEBEBE",
                                "#BEBEBE", "#BEBEBE", "#BEBEBE", "#BEBEBE", "#BEBEBE",
                                "#CC79A7", "#BEBEBE", "#009E73",
                                "#0072B2", "#BEBEBE", "#BEBEBE")) +
  scale_fill_manual(values = c("#BEBEBE", "#BEBEBE", "#BEBEBE", "#BEBEBE",
                               "#BEBEBE", "#BEBEBE", "#BEBEBE", "#BEBEBE", "#BEBEBE",
                               "#CC79A7", "#BEBEBE", "#009E73",
                               "#0072B2", "#BEBEBE", "#BEBEBE")) +
  annotate("text", label = "Lymphocyte", x = 195, y = 25, size = 8,  family = "Helvetica",
           color = "#009E73", hjust = 1) +
  annotate("text", label = "Monocyte", x = 195, y = 7.5, size = 8,  family = "Helvetica",
           color = "#CC79A7", hjust = 1) +
  annotate("text", label = "WBC", x = 195, y = 4, size = 8,  family = "Helvetica",
           color = "#0072B2", hjust = 1) +
  annotate("text", label = "12 other traits", x = 195, y = -1, size = 8,  family = "Helvetica",
           color = "#BEBEBE", hjust = 1) +
  coord_cartesian(xlim = c(upper, 200)) +
  scale_y_continuous(breaks = c(1, 10, 20, 30))

# Combine all the panels
lymphocyte_plot <- plot_grid(NULL, NULL, plot_heterozygosity, plot_var_pgs,
                            labels = c('B.', 
                                       'C.',
                                       '',
                                       ''), ncol = 2, nrow = 2,
                            label_x = 0.01, hjust = 0,
                            label_size = 28, scale = 1,
                            rel_heights = c(0.1, 1),
                            rel_widths = c(3, 2),
                            label_fontfamily = "Helvetica")

lymphocyte_plot_2 <- plot_grid(NULL, plot_compare_beta,
                              labels = c('A.', ''), ncol = 1, nrow = 2,
                              label_x = 0.01, hjust = 0,
                              label_size = 28, scale = 1,
                              rel_heights = c(0.1, 1),
                              label_fontfamily = "Helvetica")

lymphocyte_plot <- plot_grid(lymphocyte_plot_2, lymphocyte_plot,
                            nrow = 2)

grDevices::cairo_pdf("img/fig_3_lymphocyte.pdf", width = 18, height = 12)
grid.arrange(arrangeGrob(lymphocyte_plot,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

# Heterozygosity plot for supp. fig
heterozygosity <- read_tsv("data/heterozygosity/mean_heterozygosity.tsv")

`%notin%` <- Negate(`%in%`)

heterozygosity <- heterozygosity %>% 
  filter(phenotype %notin% c("Lymphocyte", "Height", "Triglycerides"))
heterozygosity$size <- ifelse(heterozygosity$size == 1, "small",
                              ifelse(heterozygosity$size == 2, "medium", "large"))
heterozygosity$size <- factor(heterozygosity$size, levels = c("small", "medium", "large"))

heterozygosity <- heterozygosity %>% filter(median_pc > upper)
heterozygosity$phenotype <- factor(heterozygosity$phenotype,
                                  levels = c("Weight", "BMI", "Body_Fat_Perc",
                                             "Monocyte", "WBC", "Eosinophil",
                                             "MCV", "MCH", "RBC",
                                             "Cystatin_C", "Platelet", "LDL"))

heterozygosity_label <- cbind.data.frame(phenotype = rep(c("Weight", "BMI", "Body fat perc",
                                                          "Monocyte", "WBC", "Eosinophil",
                                                          "MCV", "MCH", "RBC",
                                                          "Cystatin C", "Platelet", "LDL"), each = 3),
                                        size = rep(c("small", "medium", "large"), 12),
                                        label = rep(c("small effect",
                                                      "medium effect",
                                                      "large effect"), 12),
                                        median_pc = 200,
                                        mean_het = c(0.44, 0.22, 0.06, 
                                                     0.44, 0.24, 0.07,
                                                     0.44, 0.23, 0.07,
                                                     0.42, 0.16, 0.02,
                                                     0.43, 0.20, 0.04,
                                                     0.43, 0.19, 0.04,
                                                     0.41, 0.14, 0.03,
                                                     0.41, 0.14, 0.03,
                                                     0.42, 0.17, 0.04,
                                                     0.41, 0.18, 0.04,
                                                     0.42, 0.17, 0.03,
                                                     0.41, 0.14, 0.03))
heterozygosity_label$phenotype <- factor(heterozygosity_label$phenotype,
                                        levels = c("Weight", "BMI", "Body fat perc",
                                                   "Monocyte", "WBC", "Eosinophil",
                                                   "MCV", "MCH", "RBC",
                                                   "Cystatin C", "Platelet", "LDL"),
                                        labels = c("Weight", "BMI", "Body fat percentage",
                                                   "Monocyte", "White blood cell", "Eosinophil",
                                                   "Mean corpuscular volume", "Mean corpuscular hemoglobin", "Red blood cell",
                                                   "Cystatin C", "Platelet", "LDL cholesterol"))

heterozygosity_segment <- heterozygosity_label
heterozygosity_segment$median_pc <- 150
heterozygosity_segment$y <- c(0.42, 0.24, 0.08,
                              0.42, 0.26, 0.09,
                              0.42, 0.25, 0.09,
                              0.40, 0.18, 0.04,
                              0.41, 0.22, 0.06,
                              0.41, 0.21, 0.06,
                              0.39, 0.16, 0.05,
                              0.39, 0.16, 0.05,
                              0.40, 0.19, 0.06,
                              0.39, 0.20, 0.06,
                              0.40, 0.19, 0.05,
                              0.39, 0.16, 0.05)
heterozygosity_segment$yend <- c(0.39, 0.27, 0.11,
                                0.39, 0.29, 0.12,
                                0.39, 0.28, 0.12,
                                0.37, 0.21, 0.07,
                                0.38, 0.25, 0.09,
                                0.38, 0.24, 0.09,
                                0.36, 0.19, 0.08,
                                0.36, 0.19, 0.08,
                                0.37, 0.22, 0.09,
                                0.36, 0.23, 0.09,
                                0.37, 0.22, 0.08,
                                0.36, 0.19, 0.08)

plot_heterozygosity_supp <- heterozygosity %>%
  mutate(phenotype = str_replace_all(phenotype, "_", " ")) %>%
  mutate(phenotype = ifelse(phenotype == "Body Fat Perc", "Body fat perc", phenotype)) %>%
  mutate(phenotype = factor(phenotype, levels = c("Weight", "BMI", "Body fat perc",
                                                  "Monocyte", "WBC", "Eosinophil",
                                                  "MCV", "MCH", "RBC",
                                                  "Cystatin C", "Platelet", "LDL"),
                            labels = c("Weight", "BMI", "Body fat percentage",
                                       "Monocyte", "White blood cell", "Eosinophil",
                                       "Mean corpuscular volume", "Mean corpuscular hemoglobin", "Red blood cell",
                                       "Cystatin C", "Platelet", "LDL cholesterol"))) %>%
  ggplot(aes(x = median_pc, y = mean_het, color = size, fill = size)) +
  scale_x_continuous(breaks=c(0, 40, 80, 120, 160),
                     expand = c(0, 0)) +
  geom_point(size=5, alpha=0.4, shape = 23) +
  guides(color = guide_legend(title="Effect size sq.", nrow = 1),
         fill = guide_legend(title="Effect size sq.", nrow = 1)) +
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title = element_text(size=24, family = "Helvetica"),
        axis.text = element_text(size=20, family = "Helvetica"),
        strip.background = element_rect(fill="white"),
        strip.text = element_text(size=16, face = "bold", family = "Helvetica"),
        legend.position = "none") +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Heterozygosity") +
  geom_text(data = heterozygosity_label, aes(label = label), hjust = 1, size = 8,
            family = "Helvetica") +
  geom_segment(data = heterozygosity_segment, aes(xend = median_pc, y = y, yend = yend), size = 1) +
  facet_wrap(~phenotype, ncol = 3) +
  scale_color_manual(values = c("#ffb581", "#ff8934", "#e76100")) +
  scale_fill_manual(values = c("#ffb581", "#ff8934", "#e76100")) +
  coord_cartesian(xlim = c(upper, 200)) +
  guides(fill = guide_legend(title = "effect"),
         color = guide_legend(title = "effect"))

grDevices::cairo_pdf("img/fig_s17_heterozygosity.pdf", width = 18, height = 20)
print(plot_heterozygosity_supp)
dev.off()