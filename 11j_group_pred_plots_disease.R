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
pheno <- c("Alzheimer","T2D","Asthma")

# Read in group level and individual PGS files for covariates
group_pgs_df <- read_tsv("data/pgs_pred/group_pgs_df_disease.tsv") %>%
  filter(threshold == 1)
# This is not the file for diease traits
# Only using it to get the positions of knots
ind_pgs_df <- read_tsv("data/pgs_pred/ind_pgs_df.tsv") %>%
  select(-(16:34))

# Upper limit for the GWAS set
pc_dist_gwas <- read_tsv("data/pca/pc_dist_best_gwas_std.tsv")
upper <- unname(quantile(pc_dist_gwas$pc_dist, 0.975))

# Get positions of the knots by density
temp <- ind_pgs_df %>% filter(phenotype == "Height" & pc_dist > upper)
temp <- temp %>% arrange(pc_dist)
knots <- c(upper,
           temp$pc_dist[round(nrow(temp) / 9)], 
           temp$pc_dist[round(nrow(temp) / 9 * 2)], 
           temp$pc_dist[round(nrow(temp) / 9 * 3)], 
           temp$pc_dist[round(nrow(temp) / 9 * 4)],
           temp$pc_dist[round(nrow(temp) / 9 * 5)], 
           temp$pc_dist[round(nrow(temp) / 9 * 6)], 
           temp$pc_dist[round(nrow(temp) / 9 * 7)], 
           temp$pc_dist[round(nrow(temp) / 9 * 8)])

alzheimer_knots = c(upper, 
                    temp$pc_dist[round(nrow(temp) / 5)], 
                    temp$pc_dist[round(nrow(temp) / 5 * 2)], 
                    temp$pc_dist[round(nrow(temp) / 5 * 3)], 
                    temp$pc_dist[round(nrow(temp) / 5 * 4)])

# Group level plot
plot_group_level <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  # But for Alzheimer's use 5 bins as reference because it only has 25 bins
  if(trait == "Alzheimer"){
    denominator = plot_df %>%
      filter(group_close_to_gwas <= 5) %>%
      dplyr::summarise(mean_precision = mean(precision),
                       mean_recall = mean(recall))
  } else{
    denominator = plot_df %>%
      filter(group_close_to_gwas <= 50) %>%
      dplyr::summarise(mean_precision = mean(precision),
                       mean_recall = mean(recall))
  }
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  plot_df <- plot_df %>% pivot_longer(cols = precision:recall, names_to = "type", values_to = "value")
  
  # Fit a spline for both precision and recall
  temp <- subset(plot_df, type == "precision")
  if(trait == "Alzheimer"){
    lm_precision <- lm(value ~ splines::bs(median_pc, knots = alzheimer_knots), data = temp)
  } else{
    lm_precision <- lm(value ~ splines::bs(median_pc, knots = knots), data = temp)
  }
  plot_df_2 <- cbind.data.frame(median_pc = seq(1.908283, 197.5882, by = 0.4891998))
  plot_df_2$fitted_val <- unname(predict(lm_precision, newdata = plot_df_2))
  
  temp <- subset(plot_df, type == "recall")
  
  if(trait == "Alzheimer"){
    lm_recall <- lm(value ~ splines::bs(median_pc, knots = alzheimer_knots), data = temp)
  } else{
    lm_recall <- lm(value ~ splines::bs(median_pc, knots = knots), data = temp)
  }
  plot_df_3 <- cbind.data.frame(median_pc = seq(1.908283, 197.5882, by = 0.4891998))
  plot_df_3$fitted_val = unname(predict(lm_recall, newdata = plot_df_2))
  
  plot_df_2$type <- "precision"
  plot_df_3$type <- "recall"
  plot_df_2 <- plot_df_2 %>% rbind.data.frame(plot_df_3)
  
  plot <- ggplot(plot_df, aes(x = median_pc, y = value,
                              color=type, fill = type)) +
    geom_point(size = 5,alpha = 0.4, shape = 23)+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2,
              aes(x = median_pc, y = fitted_val), size = 2.5) +
    scale_fill_manual(values = c("#FFC20A", "#0C7BDC"),
                      labels = c("Precision", "Recall")) +
    scale_color_manual(values = c("#FFC20A", "#0C7BDC"),
                       labels = c("Precision", "Recall"))
  
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
    ylab("Measure") +
    scale_y_continuous(expand = c(0, 0),
                       labels = scales::number_format(accuracy = 0.1)) + 
    coord_cartesian(xlim = c(upper, 200), ylim = c(-0.04, 1.04)) +
    theme(plot.caption=element_text(size=16, family = "Helvetica"),
          legend.text=element_text(size=16, family = "Helvetica"),
          legend.position = c(0.88, 0.95),
          legend.key = element_rect(fill = NA),
          legend.title = element_blank())
  
  # Add some labels
  if(trait == "Alzheimer"){
    plot <- plot + annotate(geom = "segment", x = 39, y = 0.25, xend = 46, yend = 0.25, 
                           color = "black", size = 1) +
      annotate("text", label = "A bin of ~5172\nindividuals", x = 66, y = 0.25, size = 8, 
               family = "Helvetica")
      
  } else if(trait == "T2D"){
    plot = plot + annotate(geom = "segment", x = 105, y = 0.5, xend = 112, yend = 0.5, 
                           color = "black", size = 1) +
      annotate("text", label = "A bin of ~260\nindividuals", x = 132, y = 0.5, size = 8, 
               family = "Helvetica") 
  } else{
    plot = plot + annotate(geom = "segment", x = 107, y = 0.72, xend = 114, yend = 0.72, 
                           color = "black", size = 1) +
      annotate("text", label = "A bin of ~260\nindividuals", x = 134, y = 0.72, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Make the plots
plot_group_alzheimer <- plot_group_level(group_pgs_df, "Alzheimer", upper)

plot_group_t2d <- plot_group_level(group_pgs_df, "T2D", upper)

plot_group_asthma <- plot_group_level(group_pgs_df, "Asthma", upper)

# Combine the panels
# Main fig
fig_rev_1 <- plot_grid(NULL, NULL, plot_group_alzheimer, plot_group_t2d, 
                       NULL, NULL, plot_group_asthma, NULL, 
                       labels = c('A. Alzheimer’s disease (74 positives) (group level)', 
                                  'B. Type 2 diabetes (5,042 positives) (group level)',
                                  '',
                                  '',
                                  'C. Asthma (4,496 positives) (group level)', 
                                  '', 
                                  '',
                                  ''), ncol = 2, nrow = 4,
                       label_x = 0.01, hjust = 0,
                       label_size = 28, scale = 1,
                       rel_heights = c(0.1, 1, 0.1, 1),
                       label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_5_group_pred_disease.pdf", width = 24, height = 12)
grid.arrange(arrangeGrob(fig_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()