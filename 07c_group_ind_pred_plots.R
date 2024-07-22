library(dplyr)
library(tidyr)
library(data.table)
library(readr)
library(stringr)
library(ggplot2)
library(scales)
library(cowplot)
library(grid)
library(gridExtra)

# Phenotypes
pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
           "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
           "Body_Fat_Perc")

# Read in group level and individual PGS files for covariates
group_pgs_df <- read_tsv("data/pgs_pred/group_pgs_df.tsv") %>%
  filter(threshold == 1)
ind_pgs_df <- read_tsv("data/pgs_pred/ind_pgs_df.tsv") %>%
  select(-(16:34))

# Upper limit for the GWAS set
pc_dist_gwas <- read_tsv("data/pca/pc_dist_best_gwas_std.tsv")
upper <- unname(quantile(pc_dist_gwas$pc_dist, 0.975))

# Get positions of the knots by density
temp <- ind_pgs_df %>% filter(phenotype == "Height" & pc_dist > upper)
temp <- temp %>% arrange(pc_dist)
knots <- c(temp$pc_dist[round(nrow(temp) / 9)], 
          temp$pc_dist[round(nrow(temp) / 9 * 2)], 
          temp$pc_dist[round(nrow(temp) / 9 * 3)], 
          temp$pc_dist[round(nrow(temp) / 9 * 4)],
          temp$pc_dist[round(nrow(temp) / 9 * 5)], 
          temp$pc_dist[round(nrow(temp) / 9 * 6)], 
          temp$pc_dist[round(nrow(temp) / 9 * 7)], 
          temp$pc_dist[round(nrow(temp) / 9 * 8)])

# Group level plot
plot_group_level <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator = plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 50) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm = lm(relative_performance ~ 
            splines::bs(median_pc, knots = knots), data = plot_df)
  plot_df_2 = cbind.data.frame(median_pc = seq(1.908283, 197.5882, by = 0.4891998))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = median_pc, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2,
              aes(x = median_pc, y = fitted_val), size = 2.5, color = "#e66100")
  
  # Set the upper limit of the y-axis
  if(trait == "Height"){
    ymax <- 1.8
  }  else if(trait %in% c("Body_Fat_Perc", "Weight")){
    ymax <- 4.2
  } else if(trait == "BMI"){
    ymax <- 3.8
  } else if(trait %in% c("Eosinophil", "LDL")){
    ymax <- 3.4
  } else if(trait %in% c("MCH", "MCV")){
    ymax <- 1.8
  } else if(trait == "Platelet"){
    ymax <- 2.2
  } else if(trait %in% c("Triglycerides", "Monocyte", "Lymphocyte", "WBC", "Cystatin_C")){
    ymax <- 2.8
  } else{
    ymax <- 2.2
  }
  
  plot <- plot + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(size=24, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica")) +
    ylab(bquote("Prediction accuracy (partial"~R^2*")")) +
    scale_y_continuous(expand = c(0, 0),
                       labels = scales::number_format(accuracy = 0.1)) + 
    coord_cartesian(xlim = c(upper, 200), ylim = c(0.0, ymax))
  
  # Add some labels for height
  if(trait == "Height"){
    plot <- plot +
      annotate("text", label = "↑", x = 150, y = 1.694, size = 15, hjust = 1,
             family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 192, y = 1.694, size = 8, hjust = 1, 
               vjust = 0.8, family = "Helvetica") +
      annotate(geom = "segment", x = 108, y = 0.7, xend = 115, yend = 0.7, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of ~260\nindividuals", x = 135, y = 0.7, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
reverselog_trans <- function(base = exp(1)) {
  trans <- function(x) -log(x, base)
  inv <- function(x) base^(-x)
  trans_new(paste0("reverselog-", format(base)), trans, inv, 
            log_breaks(base = base), 
            domain = c(1e-100, Inf))
}

plot_ind_level <- function(pgs_df, trait, upper, full_range = F){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator = plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 50) %>%
    dplyr::summarise(ref_pred_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_pred_error) %>%
    select(-ref_pred_error) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(pc_dist > upper)
  
  # Fit a spline
  lm = lm(relative_performance ~ 
            splines::bs(pc_dist, knots = knots), data = plot_df)
  plot_df_2 = cbind.data.frame(pc_dist = seq(1.908283, 197.5882, by = 0.4891998))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "gray")+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2,
              aes(x = pc_dist, y = fitted_val), size=2.5, color = "#5D3A9B")
  
  plot <- plot + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(size=24, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica")) +
    ylab("Squared prediction error")
  
  # Cropped range
  if(full_range == F){
    # Set the lower and upper limits of the y-axis
    if(trait %in% c("MCV", "MCH", "Eosinophil", "Cystatin_C")){
      ymin <- 0.4
      ymax <- 4.2
    } else if(trait %in% c("Lymphocyte", "WBC")){
      ymin <- 0.18
      ymax <- 2.2
    } else{
      ymin <- 0.4
      ymax <- 2.2
    }
    
    plot <- plot +
      scale_y_continuous(expand = c(0, 0),
                         breaks = c(0.25, 0.5, 1, 2, 4),
                         labels = c("0.25", "0.50", "1.00", "2.00", "4.00"),
                         trans = reverselog_trans(2)) +
      coord_cartesian(xlim = c(upper, 200), ylim = c(ymax, ymin))
    
    # Add some labels for height
    if(trait == "Height"){
      plot <- plot +
        annotate(geom = "segment", x = 108, y = 0.75, xend = 115, yend = 0.75, 
                 color = "black", size = 1) +
        annotate("text", label = "An individual", x = 135, y = 0.75, size = 8, 
                 family = "Helvetica") +
        annotate("text", label = "↑", x = 150, y = 0.45, size = 15, hjust = 1,
                 family = "Helvetica") +
        annotate("text", label = "Better prediction", x = 192, y = 0.45, size = 8, hjust = 1, vjust = 0.8,
                 family = "Helvetica")
    }
  # Full range
  } else{
    plot <- plot +
      scale_y_continuous(expand = c(0, 0),
                         breaks=c(0.0625, 0.25, 1, 4, 16, 64, 256, 1024),
                         labels = c("0.0625", "0.25", "1", "4", "16", "64", "256", "1024"),
                         trans = reverselog_trans(2)) +
      coord_cartesian(xlim = c(upper, 200),
                      ylim = c(max(plot_df$relative_performance, na.rm = T), 0.05))
  }
  return(plot)
}

# Make the plots
plot_group_height <- plot_group_level(group_pgs_df, "Height", upper)
plot_ind_height <- plot_ind_level(ind_pgs_df, "Height", upper, full_range = F)
plot_ind_full_height <- plot_ind_level(ind_pgs_df, "Height", upper, full_range = T)

plot_group_cystatin_c <- plot_group_level(group_pgs_df, "Cystatin_C", upper)
plot_ind_cystatin_c <- plot_ind_level(ind_pgs_df, "Cystatin_C", upper, full_range = F)
plot_ind_full_cystatin_c <- plot_ind_level(ind_pgs_df, "Cystatin_C", upper, full_range = T)

plot_group_platelet <- plot_group_level(group_pgs_df, "Platelet", upper)
plot_ind_platelet <- plot_ind_level(ind_pgs_df, "Platelet", upper, full_range = F)
plot_ind_full_platelet <- plot_ind_level(ind_pgs_df, "Platelet", upper, full_range = T)

plot_group_mcv <- plot_group_level(group_pgs_df, "MCV", upper)
plot_ind_mcv <- plot_ind_level(ind_pgs_df, "MCV", upper, full_range = F)
plot_ind_full_mcv <- plot_ind_level(ind_pgs_df, "MCV", upper, full_range = T)

plot_group_weight <- plot_group_level(group_pgs_df, "Weight", upper)
plot_ind_weight <- plot_ind_level(ind_pgs_df, "Weight", upper, full_range = F)
plot_ind_full_weight <- plot_ind_level(ind_pgs_df, "Weight", upper, full_range = T)

plot_group_mch <- plot_group_level(group_pgs_df, "MCH", upper)
plot_ind_mch <- plot_ind_level(ind_pgs_df, "MCH", upper, full_range = F)
plot_ind_full_mch <- plot_ind_level(ind_pgs_df, "MCH", upper, full_range = T)

plot_group_bmi <- plot_group_level(group_pgs_df, "BMI", upper)
plot_ind_bmi <- plot_ind_level(ind_pgs_df, "BMI", upper, full_range = F)
plot_ind_full_bmi <- plot_ind_level(ind_pgs_df, "BMI", upper, full_range = T)

plot_group_rbc <- plot_group_level(group_pgs_df, "RBC", upper)
plot_ind_rbc <- plot_ind_level(ind_pgs_df, "RBC", upper, full_range = F)
plot_ind_full_rbc <- plot_ind_level(ind_pgs_df, "RBC", upper, full_range = T)

plot_group_body_fat_perc <- plot_group_level(group_pgs_df, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc <- plot_ind_level(ind_pgs_df, "Body_Fat_Perc", upper, full_range = F)
plot_ind_full_body_fat_perc <- plot_ind_level(ind_pgs_df, "Body_Fat_Perc", upper, full_range = T)

plot_group_monocyte <- plot_group_level(group_pgs_df, "Monocyte", upper)
plot_ind_monocyte <- plot_ind_level(ind_pgs_df, "Monocyte", upper, full_range = F)
plot_ind_full_monocyte <- plot_ind_level(ind_pgs_df, "Monocyte", upper, full_range = T)

plot_group_triglycerides <- plot_group_level(group_pgs_df, "Triglycerides", upper)
plot_ind_triglycerides <- plot_ind_level(ind_pgs_df, "Triglycerides", upper, full_range = F)
plot_ind_full_triglycerides <- plot_ind_level(ind_pgs_df, "Triglycerides", upper, full_range = T)

plot_group_lymphocyte <- plot_group_level(group_pgs_df, "Lymphocyte", upper)
plot_ind_lymphocyte <- plot_ind_level(ind_pgs_df, "Lymphocyte", upper, full_range = F)
plot_ind_full_lymphocyte <- plot_ind_level(ind_pgs_df, "Lymphocyte", upper, full_range = T)

plot_group_wbc <- plot_group_level(group_pgs_df, "WBC", upper)
plot_ind_wbc <- plot_ind_level(ind_pgs_df, "WBC", upper, full_range = F)
plot_ind_full_wbc <- plot_ind_level(ind_pgs_df, "WBC", upper, full_range = T)

plot_group_eosinophil <- plot_group_level(group_pgs_df, "Eosinophil", upper)
plot_ind_eosinophil <- plot_ind_level(ind_pgs_df, "Eosinophil", upper, full_range = F)
plot_ind_full_eosinophil <- plot_ind_level(ind_pgs_df, "Eosinophil", upper, full_range = T)

plot_group_ldl <- plot_group_level(group_pgs_df, "LDL", upper)
plot_ind_ldl <- plot_ind_level(ind_pgs_df, "LDL", upper, full_range = F)
plot_ind_full_ldl <- plot_ind_level(ind_pgs_df, "LDL", upper, full_range = T)

# Combine the panels
# Main fig
fig_2 <- plot_grid(NULL, NULL, plot_group_height, plot_ind_height, 
                  NULL, NULL, plot_group_weight, plot_ind_weight, 
                  NULL, NULL, plot_group_triglycerides, plot_ind_triglycerides, 
                  labels = c('A. Height (group level)', 
                             'B. Height (individual level)', 
                             '',
                             '',
                             'C. Weight (group level)', 
                             'D. Weight (individual level)', 
                             '',
                             '',
                             'E. Triglycerides (group level)', 
                             'F. Triglycerides (individual level)',
                             '',
                             ''), ncol = 2, nrow = 6,
                  label_x = 0.01, hjust = 0,
                  label_size = 28, scale = 1,
                  rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1),
                  label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_2_group_ind_pred.pdf", width = 24, height = 18)
grid.arrange(arrangeGrob(fig_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

# Supp fig for group and individual pred
fig_s2 <- plot_grid(NULL, NULL, plot_group_bmi, plot_ind_bmi, 
                    NULL, NULL, plot_group_body_fat_perc, plot_ind_body_fat_perc, 
                    labels = c('A. BMI (group level)', 
                               'B. BMI (individual level)',
                               '',
                               '',
                               'C. Body fat perc (group level)', 
                               'D. Body fat perc (individual level)', 
                               '',
                               ''), ncol = 2, nrow = 4,
                    label_x = 0.01, hjust = 0,
                    label_size = 28, scale = 1,
                    rel_heights = c(0.1, 1, 0.1, 1),
                    label_fontfamily = "Helvetica")

fig_s3 <- plot_grid(NULL, NULL, plot_group_monocyte, plot_ind_monocyte, 
                    NULL, NULL, plot_group_lymphocyte, plot_ind_lymphocyte, 
                    NULL, NULL, plot_group_wbc, plot_ind_wbc, 
                    NULL, NULL, plot_group_eosinophil, plot_ind_eosinophil, 
                    labels = c('A. Monocyte (group level)', 
                               'B. Monocyte (individual level)',
                               '',
                               '',
                               'C. Lymphocyte (group level)', 
                               'D. Lymphocyte (individual level)', 
                               '',
                               '',
                               'E. White blood cell (group level)', 
                               'F. White blood cell (individual level)', 
                               '',
                               '',
                               'G. Eosinophil (group level)', 
                               'H. Eosinophil (individual level)', 
                               '',
                               ''), ncol = 2, nrow = 8,
                    label_x = 0.01, hjust = 0,
                    label_size = 28, scale = 1,
                    rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                    label_fontfamily = "Helvetica")

fig_s4 <- plot_grid(NULL, NULL, plot_group_mcv, plot_ind_mcv, 
                    NULL, NULL, plot_group_mch, plot_ind_mch, 
                    NULL, NULL, plot_group_rbc, plot_ind_rbc, 
                    labels = c('A. Mean corpuscular volume (group level)', 
                               'B. Mean corpuscular volume (individual level)',
                               '',
                               '',
                               'C. Mean corpuscular hemoglobin (group level)', 
                               'D. Mean corpuscular hemoglobin (individual level)',
                               '',
                               '',
                               'E. Red blood cell (group level)', 
                               'F. Red blood cell (individual level)', 
                               '',
                               ''), ncol = 2, nrow = 6,
                    label_x = 0.01, hjust = 0,
                    label_size = 28, scale = 1,
                    rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1),
                    label_fontfamily = "Helvetica")

fig_s5 <- plot_grid(NULL, NULL, plot_group_cystatin_c, plot_ind_cystatin_c, 
                    NULL, NULL, plot_group_platelet, plot_ind_platelet, 
                    NULL, NULL, plot_group_ldl, plot_ind_ldl, 
                    labels = c('A. Cystatic C (group level)', 
                               'B. Cystatic C (individual level)', 
                               '',
                               '',
                               'C. Platelet (group level)', 
                               'D. Platelet (individual level)',
                               '',
                               '',
                               'E. LDL cholesterol (group level)', 
                               'F. LDL cholesterol (individual level)', 
                               '',
                               ''), ncol = 2, nrow = 6,
                    label_x = 0.01, hjust = 0,
                    label_size = 28, scale = 1,
                    rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1),
                    label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s2_group_ind_pred_physical.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_s2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s3_group_ind_pred_wbc.pdf", width = 24, height = 24, onefile = T)
grid.arrange(arrangeGrob(fig_s3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s4_group_ind_pred_rbc.pdf", width = 24, height = 18, onefile = T)
grid.arrange(arrangeGrob(fig_s4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s5_group_ind_pred_other.pdf", width = 24, height = 18, onefile = T)
grid.arrange(arrangeGrob(fig_s5,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

# Supp fig for full-ranged individual pred
fig_s6_full = plot_grid(NULL, NULL, plot_ind_full_height, plot_ind_full_weight, 
                        NULL, NULL, plot_ind_full_bmi, plot_ind_full_body_fat_perc, 
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
                        rel_heights = c(0.1, 1, 0.1, 1),
                        label_fontfamily = "Helvetica")

fig_s7_full = plot_grid(NULL, NULL, plot_ind_full_monocyte, plot_ind_full_lymphocyte, 
                        NULL, NULL, plot_ind_full_wbc, plot_ind_full_eosinophil, 
                        labels = c('A. Monocyte', 
                                   'B. Lymphocyte',
                                   '',
                                   '',
                                   'C. White blood cell', 
                                   'D. Eosinophil', 
                                   '',
                                   ''), ncol = 2, nrow = 4,
                        label_x = 0.01, hjust = 0,
                        label_size = 28, scale = 1,
                        rel_heights = c(0.1, 1, 0.1, 1),
                        label_fontfamily = "Helvetica")

fig_s8_full = plot_grid(NULL, NULL, plot_ind_full_mcv, plot_ind_full_mch, 
                        NULL, NULL, plot_ind_full_rbc, NULL, 
                        labels = c('A. Mean corpuscular volume', 
                                   'B. Mean corpuscular hemoglobin',
                                   '',
                                   '',
                                   'C. Red blood cell', 
                                   '',
                                   '',
                                   ''), ncol = 2, nrow = 4,
                        label_x = 0.01, hjust = 0,
                        label_size = 28, scale = 1,
                        rel_heights = c(0.1, 1, 0.1, 1),
                        label_fontfamily = "Helvetica")

fig_s9_full = plot_grid(NULL, NULL, plot_ind_full_cystatin_c, plot_ind_full_platelet, 
                        NULL, NULL, plot_ind_full_triglycerides, plot_ind_full_ldl, 
                        labels = c('A. Cystatin C', 
                                   'B. Platelet',
                                   '',
                                   '',
                                   'C. Triglycerides', 
                                   'D. LDL cholesterol', 
                                   '',
                                   ''), ncol = 2, nrow = 4,
                        label_x = 0.01, hjust = 0,
                        label_size = 28, scale = 1,
                        rel_heights = c(0.1, 1, 0.1, 1),
                        label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s6_ind_full_pred_physical.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_s6_full,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s7_ind_full_pred_wbc.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_s7_full,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s8_ind_full_pred_rbc.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_s8_full,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s9_ind_full_pred_other.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(fig_s9_full,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()