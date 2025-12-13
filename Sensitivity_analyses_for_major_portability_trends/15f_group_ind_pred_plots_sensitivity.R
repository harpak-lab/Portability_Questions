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
library(mgcv)

# Phenotypes
pheno <- c("Height","Platelet","MCV","MCH","BMI","RBC","Monocyte",
           "Lymphocyte","WBC","Eosinophil", "LDL", "Weight", "Triglycerides", "Cystatin_C", 
           "Body_Fat_Perc")

# Read in group level and individual PGS files for covariates
group_pgs_df <- read_tsv("data/pgs_pred/group_pgs_df.tsv") %>%
  filter(threshold == 1)
ind_pgs_df <- read_tsv("data/pgs_pred/ind_pgs_df.tsv")

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

# Original version of plots for comparison
# Group level plot
plot_group_level <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots), data = plot_df)
           
  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
  
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
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
    #ylab(bquote("Prediction accuracy (partial"~R^2*")")) +
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
      annotate("text", label = "A bin of 278\nindividuals", x = 135, y = 0.7, size = 8, 
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

plot_ind_level <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
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
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots), data = plot_df)
  plot_df_2 <- cbind.data.frame(pc_dist = seq(1.908283, 197.5882, by = 0.4891998))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2,
              aes(x = pc_dist, y = fitted_val), size=2.5, color = "#5D3A9B")
  
  plot <- plot + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(size=24, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica")) #+
    #ylab("Squared prediction error")
  
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
  return(plot)
}

# Make the plots
plot_group_height <- plot_group_level(group_pgs_df, "Height", upper)
plot_ind_height <- plot_ind_level(ind_pgs_df, "Height", upper)

plot_group_cystatin_c <- plot_group_level(group_pgs_df, "Cystatin_C", upper)
plot_ind_cystatin_c <- plot_ind_level(ind_pgs_df, "Cystatin_C", upper)

plot_group_platelet <- plot_group_level(group_pgs_df, "Platelet", upper)
plot_ind_platelet <- plot_ind_level(ind_pgs_df, "Platelet", upper)

plot_group_mcv <- plot_group_level(group_pgs_df, "MCV", upper)
plot_ind_mcv <- plot_ind_level(ind_pgs_df, "MCV", upper)

plot_group_weight <- plot_group_level(group_pgs_df, "Weight", upper)
plot_ind_weight <- plot_ind_level(ind_pgs_df, "Weight", upper)

plot_group_mch <- plot_group_level(group_pgs_df, "MCH", upper)
plot_ind_mch <- plot_ind_level(ind_pgs_df, "MCH", upper)

plot_group_bmi <- plot_group_level(group_pgs_df, "BMI", upper)
plot_ind_bmi <- plot_ind_level(ind_pgs_df, "BMI", upper)

plot_group_rbc <- plot_group_level(group_pgs_df, "RBC", upper)
plot_ind_rbc <- plot_ind_level(ind_pgs_df, "RBC", upper)

plot_group_body_fat_perc <- plot_group_level(group_pgs_df, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc <- plot_ind_level(ind_pgs_df, "Body_Fat_Perc", upper)

plot_group_monocyte <- plot_group_level(group_pgs_df, "Monocyte", upper)
plot_ind_monocyte <- plot_ind_level(ind_pgs_df, "Monocyte", upper)

plot_group_triglycerides <- plot_group_level(group_pgs_df, "Triglycerides", upper)
plot_ind_triglycerides <- plot_ind_level(ind_pgs_df, "Triglycerides", upper)

plot_group_lymphocyte <- plot_group_level(group_pgs_df, "Lymphocyte", upper)
plot_ind_lymphocyte <- plot_ind_level(ind_pgs_df, "Lymphocyte", upper)

plot_group_wbc <- plot_group_level(group_pgs_df, "WBC", upper)
plot_ind_wbc <- plot_ind_level(ind_pgs_df, "WBC", upper)

plot_group_eosinophil <- plot_group_level(group_pgs_df, "Eosinophil", upper)
plot_ind_eosinophil <- plot_ind_level(ind_pgs_df, "Eosinophil", upper)

plot_group_ldl <- plot_group_level(group_pgs_df, "LDL", upper)
plot_ind_ldl <- plot_ind_level(ind_pgs_df, "LDL", upper)


#================================= With confidence interval =================================
# Group level plot
plot_group_level_ci <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
  predictions <- predict(lm, newdata = plot_df_2, interval = "confidence")
  
  plot_df_2$fitted_val <- unname(predictions[, "fit"])
  plot_df_2$low <- unname(predictions[, "lwr"])
  plot_df_2$up <- unname(predictions[, "upr"])
  
  plot <- ggplot(plot_df, aes(x = median_pc, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_ribbon(data = plot_df_2, aes(x = median_pc, y = fitted_val, ymin = low, ymax = up), alpha = 0.2, fill = "#e66100") +
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
      annotate("text", label = "A bin of 278\nindividuals", x = 135, y = 0.7, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_ci <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
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
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  predictions <- predict(lm, newdata = plot_df_2, interval = "confidence")
  
  plot_df_2$fitted_val <- unname(predictions[, "fit"])
  plot_df_2$low <- unname(predictions[, "lwr"])
  plot_df_2$up <- unname(predictions[, "upr"])
  
  # Replace negative or zero values with a small positive value for log2 scale
  if(any(plot_df_2$low <= 0)){
    print(trait)
  }
  plot_df_2$low[plot_df_2$low <= 0] <- 1e-6
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_ribbon(data = plot_df_2, aes(x = pc_dist, y = fitted_val, ymin = low, ymax = up), alpha = 0.2, fill = "#5D3A9B") +
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
  
  if(trait %in% c("MCV", "MCH", "Eosinophil")){
    ymin <- 0.4
    ymax <- 4.2
  } else if(trait == "WBC"){
    ymin <- 0.18
    ymax <- 2.2
  } else if(trait == "Cystatin_C"){
    ymin <- 0.02
    ymax <- 4.2
  } else if(trait %in% c("Lymphocyte", "Monocyte")){
    ymin <- 0.02
    ymax <- 2.2
  } else{
    ymin <- 0.4
    ymax <- 2.2
  }
  
  plot <- plot +
    scale_y_continuous(expand = c(0, 0),
                       breaks=c(0.03125, 0.0625, 0.125, 0.25, 0.5, 1, 2, 4),
                       labels = c("0.03125", "0.0625", "0.125", "0.25", "0.50", "1.00", "2.00", "4.00"),
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
  return(plot)
}

# Make the plots
plot_group_height_ci <- plot_group_level_ci(group_pgs_df, "Height", upper)
plot_ind_height_ci <- plot_ind_level_ci(ind_pgs_df, "Height", upper)

plot_group_cystatin_c_ci <- plot_group_level_ci(group_pgs_df, "Cystatin_C", upper)
plot_ind_cystatin_c_ci <- plot_ind_level_ci(ind_pgs_df, "Cystatin_C", upper)

plot_group_platelet_ci <- plot_group_level_ci(group_pgs_df, "Platelet", upper)
plot_ind_platelet_ci <- plot_ind_level_ci(ind_pgs_df, "Platelet", upper)

plot_group_mcv_ci <- plot_group_level_ci(group_pgs_df, "MCV", upper)
plot_ind_mcv_ci <- plot_ind_level_ci(ind_pgs_df, "MCV", upper)

plot_group_weight_ci <- plot_group_level_ci(group_pgs_df, "Weight", upper)
plot_ind_weight_ci <- plot_ind_level_ci(ind_pgs_df, "Weight", upper)

plot_group_mch_ci <- plot_group_level_ci(group_pgs_df, "MCH", upper)
plot_ind_mch_ci <- plot_ind_level_ci(ind_pgs_df, "MCH", upper)

plot_group_bmi_ci <- plot_group_level_ci(group_pgs_df, "BMI", upper)
plot_ind_bmi_ci <- plot_ind_level_ci(ind_pgs_df, "BMI", upper)

plot_group_rbc_ci <- plot_group_level_ci(group_pgs_df, "RBC", upper)
plot_ind_rbc_ci <- plot_ind_level_ci(ind_pgs_df, "RBC", upper)

plot_group_body_fat_perc_ci <- plot_group_level_ci(group_pgs_df, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_ci <- plot_ind_level_ci(ind_pgs_df, "Body_Fat_Perc", upper)

plot_group_monocyte_ci <- plot_group_level_ci(group_pgs_df, "Monocyte", upper)
plot_ind_monocyte_ci <- plot_ind_level_ci(ind_pgs_df, "Monocyte", upper)

plot_group_triglycerides_ci <- plot_group_level_ci(group_pgs_df, "Triglycerides", upper)
plot_ind_triglycerides_ci <- plot_ind_level_ci(ind_pgs_df, "Triglycerides", upper)

plot_group_lymphocyte_ci <- plot_group_level_ci(group_pgs_df, "Lymphocyte", upper)
plot_ind_lymphocyte_ci <- plot_ind_level_ci(ind_pgs_df, "Lymphocyte", upper)

plot_group_wbc_ci <- plot_group_level_ci(group_pgs_df, "WBC", upper)
plot_ind_wbc_ci <- plot_ind_level_ci(ind_pgs_df, "WBC", upper)

plot_group_eosinophil_ci <- plot_group_level_ci(group_pgs_df, "Eosinophil", upper)
plot_ind_eosinophil_ci <- plot_ind_level_ci(ind_pgs_df, "Eosinophil", upper)

plot_group_ldl_ci <- plot_group_level_ci(group_pgs_df, "LDL", upper)
plot_ind_ldl_ci <- plot_ind_level_ci(ind_pgs_df, "LDL", upper)

plot_ci_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                       plot_group_height_ci, plot_group_height, plot_ind_height_ci, plot_ind_height, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_weight_ci, plot_group_weight, plot_ind_weight_ci, plot_ind_weight,
                       NULL, NULL, NULL, NULL, 
                       plot_group_bmi_ci, plot_group_bmi, plot_ind_bmi_ci, plot_ind_bmi, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_body_fat_perc_ci, plot_group_body_fat_perc, plot_ind_body_fat_perc_ci, plot_ind_body_fat_perc, 
                       labels = c('A. Height (group level, w/ CI)',
                                  'B. Height (group level, main text ver.)',
                                  'C. Height (individual level, w/ CI)',
                                  'D. Height (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'E. Weight (group level, w/ CI)',
                                  'F. Weight (group level, main text ver.)',
                                  'G. Weight (individual level, w/ CI)',
                                  'H. Weight (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'I. BMI (group level, w/ CI)', 
                                  'J. BMI (group level, main text ver.)',
                                  'K. BMI (individual level, w/ CI)',
                                  'L. BMI (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'M. Body fat percentage (group level, w/ CI)', 
                                  'N. Body fat percentage (group level, main text ver.)', 
                                  'O. Body fat percentage (individual level, w/ CI)',
                                  'P. Body fat percentage (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                 ''), 
                       ncol = 4, nrow = 8,
                       label_x = 0.01, hjust = 0,
                       label_size = 28, scale = 1,
                       rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                       label_fontfamily = "Helvetica")

plot_ci_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                       plot_group_monocyte_ci, plot_group_monocyte, plot_ind_monocyte_ci, plot_ind_monocyte, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_lymphocyte_ci, plot_group_lymphocyte, plot_ind_lymphocyte_ci, plot_ind_lymphocyte, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_wbc_ci, plot_group_wbc, plot_ind_wbc_ci, plot_ind_wbc, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_eosinophil_ci, plot_group_eosinophil, plot_ind_eosinophil_ci, plot_ind_eosinophil, 
                       labels = c('A. Monocyte count (group level, w/ CI)', 
                                  'B. Monocyte count (group level, main text ver.)',
                                  'C. Monocyte count (individual level, w/ CI)',
                                  'D. Monocyte count (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'E. Lymphocyte count (group level, w/ CI)', 
                                  'F. Lymphocyte count (group level, main text ver.)',
                                  'G. Lymphocyte count (individual level, w/ CI)', 
                                  'H. Lymphocyte count (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'I. White blood cell count (group level, w/ CI)', 
                                  'J. White blood cell count (group level, main text ver.)',
                                  'K. White blood cell count (individual level, w/ CI)', 
                                  'L. White blood cell count (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'M. Eosinophil count (group level, w/ CI)', 
                                  'N. Eosinophil count (group level, main text ver.)',
                                  'O. Eosinophil count (individual level, w/ CI)', 
                                  'P. Eosinophil count (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  ''), 
                       ncol = 4, nrow = 8,
                       label_x = 0.01, hjust = 0,
                       label_size = 28, scale = 1,
                       rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                       label_fontfamily = "Helvetica")

plot_ci_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                       plot_group_mcv_ci, plot_group_mcv, plot_ind_mcv_ci, plot_ind_mcv, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_mch_ci, plot_group_mch, plot_ind_mch_ci, plot_ind_mch,
                       NULL, NULL, NULL, NULL, 
                       plot_group_rbc_ci, plot_group_rbc, plot_ind_rbc_ci, plot_ind_rbc, 
                       labels = c('A. Mean corpuscular volume (group level, w/ CI)', 
                                  'B. Mean corpuscular volume (group level, main text ver.)',
                                  'C. Mean corpuscular volume (individual level, w/ CI)',
                                  'D. Mean corpuscular volume (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'E. Mean corpuscular hemoglobin (group level, w/ CI)\n',
                                  'F. Mean corpuscular hemoglobin (group level, main text ver.)\n',
                                  'G. Mean corpuscular hemoglobin (individual level, w/ CI)\n',
                                  'H. Mean corpuscular hemoglobin (individual level, main text\nver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'I. Red blood cell count (group level, w/ CI)', 
                                  'J. Red blood cell count (group level, main text ver.)',
                                  'K. Red blood cell count (individual level, w/ CI)', 
                                  'L. Red blood cell count (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  ''), 
                       ncol = 4, nrow = 6,
                       label_x = 0.01, hjust = 0,
                       label_size = 28, scale = 1,
                       rel_heights = c(0.1, 1, 0.25, 1, 0.1, 1),
                       label_fontfamily = "Helvetica")

plot_ci_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                       plot_group_cystatin_c_ci, plot_group_cystatin_c, plot_ind_cystatin_c_ci, plot_ind_cystatin_c, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_platelet_ci, plot_group_platelet, plot_ind_platelet_ci, plot_ind_platelet, 
                       NULL, NULL, NULL, NULL, 
                       plot_group_triglycerides_ci, plot_group_triglycerides, plot_ind_triglycerides_ci, plot_ind_triglycerides,
                       NULL, NULL, NULL, NULL, 
                       plot_group_ldl_ci, plot_group_ldl, plot_ind_ldl_ci, plot_ind_ldl, 
                       labels = c('A. Cystatin C level (group level, w/ CI)', 
                                  'B. Cystatin C level (group level, main text ver.)',
                                  'C. Cystatin C level (individual level, w/ CI)', 
                                  'D. Cystatin C level (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'E. Platelet count (group level, w/ CI)', 
                                  'F. Platelet count (group level, main text ver.)', 
                                  'G. Platelet count (individual level, w/ CI)',
                                  'H. Platelet count (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'I. Triglyceride level (group level, w/ CI)', 
                                  'J. Triglyceride level (group level, main text ver.)', 
                                  'K. Triglyceride level (individual level, w/ CI)', 
                                  'L. Triglyceride level (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  '',
                                  'M. LDL cholesterol level (group level, w/ CI)', 
                                  'N. LDL cholesterol level (group level, main text ver.)',
                                  'O. LDL cholesterol level (individual level, w/ CI)', 
                                  'P. LDL cholesterol level (individual level, main text ver.)',
                                  '',
                                  '',
                                  '',
                                  ''), 
                       ncol = 4, nrow = 8,
                       label_x = 0.01, hjust = 0,
                       label_size = 28, scale = 1,
                       rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                       label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s34_pred_w_ci_physical.pdf", width = 48, height = 24, onefile = T)
grid.arrange(arrangeGrob(plot_ci_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s35_pred_w_ci_wbc.pdf", width = 48, height = 24, onefile = T)
grid.arrange(arrangeGrob(plot_ci_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s36_pred_w_ci_rbc.pdf", width = 48, height = 18.82, onefile = T)
grid.arrange(arrangeGrob(plot_ci_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s37_pred_w_ci_other.pdf", width = 48, height = 24, onefile = T)
grid.arrange(arrangeGrob(plot_ci_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

#================================= 500 bins =================================

# Read in group level PGS files for covariates
group_pgs_df_500_bins <- read_tsv("data/pgs_pred/group_pgs_df_500_bins.tsv") %>%
  filter(threshold == 1)

# Group level plot
plot_group_level_500_bins <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
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
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
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
  }  else if(trait %in% c("Body_Fat_Perc")){
    ymax <- 5.6
  } else if(trait %in% c("BMI", "Eosinophil")){
    ymax <- 3.8
  } else if(trait == "LDL"){
    ymax <- 3.4
  } else if(trait %in% c("MCH", "MCV")){
    ymax <- 1.8
  } else if(trait == "Platelet"){
    ymax <- 2.2
  } else if(trait == "RBC"){
    ymax <- 2.6
  } else if(trait %in% c("Triglycerides", "WBC", "Cystatin_C")){
    ymax <- 2.8
  } else if(trait %in% c("Monocyte", "Lymphocyte")){
    ymax <- 3.2
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
      annotate(geom = "segment", x = 133, y = 0.75, xend = 140, yend = 0.75, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 139\nindividuals", x = 160, y = 0.75, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Make the plots
plot_group_height_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Height", upper)

plot_group_cystatin_c_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Cystatin_C", upper)

plot_group_platelet_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Platelet", upper)

plot_group_mcv_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "MCV", upper)

plot_group_weight_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Weight", upper)

plot_group_mch_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "MCH", upper)

plot_group_bmi_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "BMI", upper)

plot_group_rbc_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "RBC", upper)

plot_group_body_fat_perc_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Body_Fat_Perc", upper)

plot_group_monocyte_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Monocyte", upper)

plot_group_triglycerides_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Triglycerides", upper)

plot_group_lymphocyte_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Lymphocyte", upper)

plot_group_wbc_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "WBC", upper)

plot_group_eosinophil_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "Eosinophil", upper)

plot_group_ldl_500_bins <- plot_group_level_500_bins(group_pgs_df_500_bins, "LDL", upper)

plot_500_bins_1 <- plot_grid(NULL, NULL, NULL, NULL,
                              plot_group_height_500_bins, plot_group_height, plot_group_weight_500_bins, plot_group_weight,
                              NULL, NULL, NULL, NULL, 
                              plot_group_bmi_500_bins, plot_group_bmi, plot_group_body_fat_perc_500_bins, plot_group_body_fat_perc,
                              labels = c('A. Height (group level, 500 bins)',
                                         'B. Height (group level, 250 bins, main text ver.)',
                                         'C. Weight (group level, 500 bins)',
                                         'D. Weight (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         '',
                                         'E. BMI (group level, 500 bins)',
                                         'F. BMI (group level, 250 bins, main text ver.)',
                                         'G. Body fat percentage (group level, 500 bins)',
                                         'H. Body fat percentage (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         ''), 
                              ncol = 4, nrow = 4,
                              label_x = 0.01, hjust = 0,
                              label_size = 28, scale = 1,
                              rel_heights = c(0.1, 1, 0.1, 1),
                              label_fontfamily = "Helvetica")

plot_500_bins_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                              plot_group_monocyte_500_bins, plot_group_monocyte, plot_group_lymphocyte_500_bins, plot_group_lymphocyte,
                              NULL, NULL, NULL, NULL, 
                              plot_group_wbc_500_bins, plot_group_wbc, plot_group_eosinophil_500_bins, plot_group_eosinophil, 
                              labels = c('A. Monocyte count (group level, 500 bins)', 
                                         'B. Monocyte count (group level, 250 bins, main text ver.)',
                                         'C. Lymphocyte count (group level, 500 bins)',
                                         'D. Lymphocyte count (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         '',
                                         'E. White blood cell count (group level, 500 bins)', 
                                         'F. White blood cell count (group level, 250 bins, main text ver.)', 
                                         'G. Eosinophil count (group level, 500 bins)',
                                         'H. Eosinophil count (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         ''), 
                              ncol = 4, nrow = 4,
                              label_x = 0.01, hjust = 0,
                              label_size = 28, scale = 1,
                              rel_heights = c(0.1, 1, 0.1, 1),
                              label_fontfamily = "Helvetica")

plot_500_bins_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                              plot_group_mch_500_bins, plot_group_mch, plot_group_mcv_500_bins, plot_group_mch, 
                              NULL, NULL, NULL, NULL, 
                              plot_group_rbc_500_bins, plot_group_rbc, NULL, NULL, 
                              labels = c('A. Mean corpuscular volume (group level, 500 bins)\n', 
                                         'B. Mean corpuscular volume (group level, 250 bins, main text\nver.)',
                                         'C. Mean corpuscular hemoglobin (group level, 500 bins)\n',
                                         'D. Mean corpuscular hemoglobin (group level, 250 bins, main\ntext ver.)',
                                         '',
                                         '',
                                         '',
                                         '',
                                         'E. Red blood cell count (group level, 500 bins)', 
                                         'F. Red blood cell count (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         ''), 
                              ncol = 4, nrow = 4,
                              label_x = 0.01, hjust = 0,
                              label_size = 28, scale = 1,
                              rel_heights = c(0.25, 1, 0.1, 1),
                              label_fontfamily = "Helvetica")

plot_500_bins_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                              plot_group_cystatin_c_500_bins, plot_group_cystatin_c, plot_group_platelet_500_bins, plot_group_platelet,
                              NULL, NULL, NULL, NULL, 
                              plot_group_triglycerides_500_bins, plot_group_triglycerides, plot_group_ldl_500_bins, plot_group_ldl, 
                              labels = c('A. Cystatin C level (group level, 500 bins)', 
                                         'B. Cystatin C level (group level, 250 bins, main text ver.)', 
                                         'C. Platelet count (group level, 500 bins)',
                                         'D. Platelet count (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         '',
                                         'E. Triglyceride level (group level, 500 bins)', 
                                         'F. Triglyceride level (group level, 250 bins, main text ver.)',
                                         'G. LDL cholesterol level (group level, 500 bins)',
                                         'H. LDL cholesterol level (group level, 250 bins, main text ver.)',
                                         '',
                                         '',
                                         '',
                                         ''), 
                              ncol = 4, nrow = 4,
                              label_x = 0.01, hjust = 0,
                              label_size = 28, scale = 1,
                              rel_heights = c(0.1, 1, 0.1, 1),
                              label_fontfamily = "Helvetica")


grDevices::cairo_pdf("img/fig_s38_500_bins_physical.pdf", width = 48, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_500_bins_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s39_500_bins_wbc.pdf", width = 48, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_500_bins_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s40_500_bins_rbc.pdf", width = 48, height = 12.81, onefile = T)
grid.arrange(arrangeGrob(plot_500_bins_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s41_500_bins_other.pdf", width = 48, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_500_bins_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= GAM =================================
# Group level plot
plot_group_level_gam <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a GAM model
  gam_model = gam(relative_performance ~ s(median_pc), knots = list(knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 = cbind.data.frame(median_pc = seq(start, stop, by = step))
  plot_df_2$fitted_val <- predict(gam_model, newdata = plot_df_2)
  
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
      annotate("text", label = "A bin of 278\nindividuals", x = 135, y = 0.7, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_gam <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_pred_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_pred_error) %>%
    select(-ref_pred_error) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(pc_dist > upper)
  
  # Fit a GAM model
  gam_model = gam(relative_performance ~ s(pc_dist), knots = list(knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 = cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val <- predict(gam_model, newdata = plot_df_2)
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
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
  return(plot)
}

# Make the plots
plot_group_height_gam <- plot_group_level_gam(group_pgs_df, "Height", upper)
plot_ind_height_gam <- plot_ind_level_gam(ind_pgs_df, "Height", upper)

plot_group_cystatin_c_gam <- plot_group_level_gam(group_pgs_df, "Cystatin_C", upper)
plot_ind_cystatin_c_gam <- plot_ind_level_gam(ind_pgs_df, "Cystatin_C", upper)

plot_group_platelet_gam <- plot_group_level_gam(group_pgs_df, "Platelet", upper)
plot_ind_platelet_gam <- plot_ind_level_gam(ind_pgs_df, "Platelet", upper)

plot_group_mcv_gam <- plot_group_level_gam(group_pgs_df, "MCV", upper)
plot_ind_mcv_gam <- plot_ind_level_gam(ind_pgs_df, "MCV", upper)

plot_group_weight_gam <- plot_group_level_gam(group_pgs_df, "Weight", upper)
plot_ind_weight_gam <- plot_ind_level_gam(ind_pgs_df, "Weight", upper)

plot_group_mch_gam <- plot_group_level_gam(group_pgs_df, "MCH", upper)
plot_ind_mch_gam <- plot_ind_level_gam(ind_pgs_df, "MCH", upper)

plot_group_bmi_gam <- plot_group_level_gam(group_pgs_df, "BMI", upper)
plot_ind_bmi_gam <- plot_ind_level_gam(ind_pgs_df, "BMI", upper)

plot_group_rbc_gam <- plot_group_level_gam(group_pgs_df, "RBC", upper)
plot_ind_rbc_gam <- plot_ind_level_gam(ind_pgs_df, "RBC", upper)

plot_group_body_fat_perc_gam <- plot_group_level_gam(group_pgs_df, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_gam <- plot_ind_level_gam(ind_pgs_df, "Body_Fat_Perc", upper)

plot_group_monocyte_gam <- plot_group_level_gam(group_pgs_df, "Monocyte", upper)
plot_ind_monocyte_gam <- plot_ind_level_gam(ind_pgs_df, "Monocyte", upper)

plot_group_triglycerides_gam <- plot_group_level_gam(group_pgs_df, "Triglycerides", upper)
plot_ind_triglycerides_gam <- plot_ind_level_gam(ind_pgs_df, "Triglycerides", upper)

plot_group_lymphocyte_gam <- plot_group_level_gam(group_pgs_df, "Lymphocyte", upper)
plot_ind_lymphocyte_gam <- plot_ind_level_gam(ind_pgs_df, "Lymphocyte", upper)

plot_group_wbc_gam <- plot_group_level_gam(group_pgs_df, "WBC", upper)
plot_ind_wbc_gam <- plot_ind_level_gam(ind_pgs_df, "WBC", upper)

plot_group_eosinophil_gam <- plot_group_level_gam(group_pgs_df, "Eosinophil", upper)
plot_ind_eosinophil_gam <- plot_ind_level_gam(ind_pgs_df, "Eosinophil", upper)

plot_group_ldl_gam <- plot_group_level_gam(group_pgs_df, "LDL", upper)
plot_ind_ldl_gam <- plot_ind_level_gam(ind_pgs_df, "LDL", upper)

plot_gam_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                        plot_group_height_gam, plot_group_height, plot_ind_height_gam, plot_ind_height, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_weight_gam, plot_group_weight, plot_ind_weight_gam, plot_ind_weight,
                        NULL, NULL, NULL, NULL, 
                        plot_group_bmi_gam, plot_group_bmi, plot_ind_bmi_gam, plot_ind_bmi, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_body_fat_perc_gam, plot_group_body_fat_perc, plot_ind_body_fat_perc_gam, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, GAM)',
                              'B. Height (group level, cubic spline, main text ver.)',
                              'C. Height (individual level, GAM)',
                              'D. Height (individual level, cubic spline, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, GAM)',
                              'F. Weight (group level, cubic spline, main text ver.)',
                              'G. Weight (individual level, GAM)',
                              'H. Weight (individual level, cubic spline, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, GAM)', 
                              'J. BMI (group level, cubic spline, main text ver.)',
                              'K. BMI (individual level, GAM)',
                              'L. BMI (individual level, cubic spline, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, GAM)\n', 
                              'N. Body fat percentage (group level, cubic spline, main text\nver.)', 
                              'O. Body fat percentage (individual level, GAM)\n',
                              'P. Body fat percentage (individual level, cubic spline, main text\nver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_gam_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                        plot_group_monocyte_gam, plot_group_monocyte, plot_ind_monocyte_gam, plot_ind_monocyte, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_lymphocyte_gam, plot_group_lymphocyte, plot_ind_lymphocyte_gam, plot_ind_lymphocyte, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_wbc_gam, plot_group_wbc, plot_ind_wbc_gam, plot_ind_wbc, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_eosinophil_gam, plot_group_eosinophil, plot_ind_eosinophil_gam, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, GAM)', 
                              'B. Monocyte count (group level, cubic spline, main text ver.)',
                              'C. Monocyte count (individual level, GAM)',
                              'D. Monocyte count (individual level, cubic spline, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, GAM)\n', 
                              'F. Lymphocyte count (group level, cubic spline, main text ver.)\n',
                              'G. Lymphocyte count (individual level, GAM)\n', 
                              'H. Lymphocyte count (individual level, cubic spline, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, GAM)\n', 
                              'J. White blood cell count (group level, cubic spline, main text\nver.)',
                              'K. White blood cell count (individual level, GAM)\n', 
                              'L. White blood cell count (individual level, cubic spline, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, GAM)\n', 
                              'N. Eosinophil count (group level, cubic spline, main text ver.)\n',
                              'O. Eosinophil count (individual level, GAM)\n', 
                              'P. Eosinophil count (individual level, cubic spline, main text\nver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_gam_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                        plot_group_mcv_gam, plot_group_mcv, plot_ind_mcv_gam, plot_ind_mcv, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_mch_gam, plot_group_mch, plot_ind_mch_gam, plot_ind_mch,
                        NULL, NULL, NULL, NULL, 
                        plot_group_rbc_gam, plot_group_rbc, plot_ind_rbc_gam, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level, GAM)\n', 
                              'B. Mean corpuscular volume (group level, cubic spline, main\ntext ver.)',
                              'C. Mean corpuscular volume (individual level, GAM)\n',
                              'D. Mean corpuscular volume (individual level, cubic spline,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, GAM)\n',
                              'F. Mean corpuscular hemoglobin (group level, cubic spline,\nmain text ver.)',
                              'G. Mean corpuscular hemoglobin (individual level, GAM)\n',
                              'H. Mean corpuscular hemoglobin (individual level, cubic spline,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, GAM)\n', 
                              'J. Red blood cell count (group level, cubic spline, main text ver.)\n',
                              'K. Red blood cell count (individual level, GAM)\n', 
                              'L. Red blood cell count (individual level, cubic spline, main text\nver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_gam_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                        plot_group_cystatin_c_gam, plot_group_cystatin_c, plot_ind_cystatin_c_gam, plot_ind_cystatin_c, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_platelet_gam, plot_group_platelet, plot_ind_platelet_gam, plot_ind_platelet, 
                        NULL, NULL, NULL, NULL, 
                        plot_group_triglycerides_gam, plot_group_triglycerides, plot_ind_triglycerides_gam, plot_ind_triglycerides,
                        NULL, NULL, NULL, NULL, 
                        plot_group_ldl_gam, plot_group_ldl, plot_ind_ldl_gam, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, GAM)', 
                              'B. Cystatin C level (group level, cubic spline, main text ver.)',
                              'C. Cystatin C level (individual level, GAM)', 
                              'D. Cystatin C level (individual level, cubic spline, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, GAM)', 
                              'F. Platelet count (group level, cubic spline, main text ver.)', 
                              'G. Platelet count (individual level, GAM)',
                              'H. Platelet count (individual level, cubic spline, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, GAM)\n', 
                              'J. Triglyceride level (group level, cubic spline, main text ver.)\n', 
                              'K. Triglyceride level (individual level, GAM)\n', 
                              'L. Triglyceride level (individual level, cubic spline, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, GAM)\n', 
                              'N. LDL cholesterol level (group level, cubic spline, main text\nver.)',
                              'O. LDL cholesterol level (individual level, GAM)\n', 
                              'P. LDL cholesterol level (individual level, cubic spline, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s42_gam_physical.pdf", width = 48, height = 24.82, onefile = T)
grid.arrange(arrangeGrob(plot_gam_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s43_gam_wbc.pdf", width = 48, height = 26.45, onefile = T)
grid.arrange(arrangeGrob(plot_gam_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s44_gam_rbc.pdf", width = 48, height = 20.45, onefile = T)
grid.arrange(arrangeGrob(plot_gam_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s45_gam_other.pdf", width = 48, height = 25.64, onefile = T)
grid.arrange(arrangeGrob(plot_gam_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= 12 knots =================================
# Get positions of the knots by density
knots_12 <- c(upper,
              temp$pc_dist[round(nrow(temp) / 13)], 
              temp$pc_dist[round(nrow(temp) / 13 * 2)], 
              temp$pc_dist[round(nrow(temp) / 13 * 3)], 
              temp$pc_dist[round(nrow(temp) / 13 * 4)],
              temp$pc_dist[round(nrow(temp) / 13 * 5)], 
              temp$pc_dist[round(nrow(temp) / 13 * 6)], 
              temp$pc_dist[round(nrow(temp) / 13 * 7)], 
              temp$pc_dist[round(nrow(temp) / 13 * 8)],
              temp$pc_dist[round(nrow(temp) / 13 * 9)],
              temp$pc_dist[round(nrow(temp) / 13 * 10)],
              temp$pc_dist[round(nrow(temp) / 13 * 11)],
              temp$pc_dist[round(nrow(temp) / 13 * 12)])

# Group level plot
plot_group_level_12_knots <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots_12), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
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
      annotate("text", label = "A bin of 278\nindividuals", x = 135, y = 0.7, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_12_knots <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
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
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots_12), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
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
  
  if(trait %in% c("MCV", "MCH", "Eosinophil")){
    ymin <- 0.4
    ymax <- 4.2
  } else if(trait %in% c("Lymphocyte", "WBC")){
    ymin <- 0.18
    ymax <- 2.2
  } else if(trait == "Lymphocyte"){
    ymin <- 0.08
    ymax <- 2.2
  } else if(trait == "Monocyte"){
    ymin <- 0.35
    ymax <- 2.2
  } else if(trait == "Cystatin_C"){
    ymin <- 0.4
    ymax <- 6.2
  } else{
    ymin <- 0.4
    ymax <- 2.2
  }
  
  plot <- plot +
    scale_y_continuous(expand = c(0, 0),
                       breaks=c(0.125, 0.25, 0.5, 1, 2, 4),
                       labels = c("0.125", "0.25", "0.50", "1.00", "2.00", "4.00"),
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
  return(plot)
}

# Make the plots
plot_group_height_12_knots <- plot_group_level_12_knots(group_pgs_df, "Height", upper)
plot_ind_height_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Height", upper)

plot_group_cystatin_c_12_knots <- plot_group_level_12_knots(group_pgs_df, "Cystatin_C", upper)
plot_ind_cystatin_c_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Cystatin_C", upper)

plot_group_platelet_12_knots <- plot_group_level_12_knots(group_pgs_df, "Platelet", upper)
plot_ind_platelet_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Platelet", upper)

plot_group_mcv_12_knots <- plot_group_level_12_knots(group_pgs_df, "MCV", upper)
plot_ind_mcv_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "MCV", upper)

plot_group_weight_12_knots <- plot_group_level_12_knots(group_pgs_df, "Weight", upper)
plot_ind_weight_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Weight", upper)

plot_group_mch_12_knots <- plot_group_level_12_knots(group_pgs_df, "MCH", upper)
plot_ind_mch_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "MCH", upper)

plot_group_bmi_12_knots <- plot_group_level_12_knots(group_pgs_df, "BMI", upper)
plot_ind_bmi_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "BMI", upper)

plot_group_rbc_12_knots <- plot_group_level_12_knots(group_pgs_df, "RBC", upper)
plot_ind_rbc_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "RBC", upper)

plot_group_body_fat_perc_12_knots <- plot_group_level_12_knots(group_pgs_df, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Body_Fat_Perc", upper)

plot_group_monocyte_12_knots <- plot_group_level_12_knots(group_pgs_df, "Monocyte", upper)
plot_ind_monocyte_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Monocyte", upper)

plot_group_triglycerides_12_knots <- plot_group_level_12_knots(group_pgs_df, "Triglycerides", upper)
plot_ind_triglycerides_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Triglycerides", upper)

plot_group_lymphocyte_12_knots <- plot_group_level_12_knots(group_pgs_df, "Lymphocyte", upper)
plot_ind_lymphocyte_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Lymphocyte", upper)

plot_group_wbc_12_knots <- plot_group_level_12_knots(group_pgs_df, "WBC", upper)
plot_ind_wbc_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "WBC", upper)

plot_group_eosinophil_12_knots <- plot_group_level_12_knots(group_pgs_df, "Eosinophil", upper)
plot_ind_eosinophil_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "Eosinophil", upper)

plot_group_ldl_12_knots <- plot_group_level_12_knots(group_pgs_df, "LDL", upper)
plot_ind_ldl_12_knots <- plot_ind_level_12_knots(ind_pgs_df, "LDL", upper)

plot_12_knots_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                             plot_group_height_12_knots, plot_group_height, plot_ind_height_12_knots, plot_ind_height, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_weight_12_knots, plot_group_weight, plot_ind_weight_12_knots, plot_ind_weight,
                             NULL, NULL, NULL, NULL, 
                             plot_group_bmi_12_knots, plot_group_bmi, plot_ind_bmi_12_knots, plot_ind_bmi, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_body_fat_perc_12_knots, plot_group_body_fat_perc, plot_ind_body_fat_perc_12_knots, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, spline w/ 12 knots)',
                              'B. Height (group level, spline w/ 8 knots, main text ver.)',
                              'C. Height (individual level, spline w/ 12 knots)',
                              'D. Height (individual level, spline w/ 8 knots, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, spline w/ 12 knots)',
                              'F. Weight (group level, spline w/ 8 knots, main text ver.)',
                              'G. Weight (individual level, spline w/ 12 knots)',
                              'H. Weight (individual level, spline w/ 8 knots, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, spline w/ 12 knots)', 
                              'J. BMI (group level, spline w/ 8 knots, main text ver.)',
                              'K. BMI (individual level, spline w/ 12 knots)',
                              'L. BMI (individual level, spline w/ 8 knots, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, spline w/ 12 knots)\n', 
                              'N. Body fat percentage (group level, spline w/ 8 knots, main text\nver.)', 
                              'O. Body fat percentage (individual level, spline w/ 12 knots)\n',
                              'P. Body fat percentage (individual level, spline w/ 8 knots, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_12_knots_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                             plot_group_monocyte_12_knots, plot_group_monocyte, plot_ind_monocyte_12_knots, plot_ind_monocyte, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_lymphocyte_12_knots, plot_group_lymphocyte, plot_ind_lymphocyte_12_knots, plot_ind_lymphocyte, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_wbc_12_knots, plot_group_wbc, plot_ind_wbc_12_knots, plot_ind_wbc, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_eosinophil_12_knots, plot_group_eosinophil, plot_ind_eosinophil_12_knots, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, spline w/ 12 knots)\n', 
                              'B. Monocyte count (group level, spline w/ 8 knots, main text\nver.)',
                              'C. Monocyte count (individual level, spline w/ 12 knots)\n',
                              'D. Monocyte count (individual level, spline w/ 8 knots, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, spline w/ 12 knots)\n', 
                              'F. Lymphocyte count (group level, spline w/ 8 knots, main text\nver.)',
                              'G. Lymphocyte count (individual level, spline w/ 12 knots)\n', 
                              'H. Lymphocyte count (individual level, spline w/ 8 knots, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, spline w/ 12 knots)\n', 
                              'J. White blood cell count (group level, spline w/ 8 knots, main\ntext ver.)',
                              'K. White blood cell count (individual level, spline w/ 12 knots)\n', 
                              'L. White blood cell count (individual level, spline w/ 8 knots,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, spline w/ 12 knots)\n', 
                              'N. Eosinophil count (group level, spline w/ 8 knots, main text\nver.)',
                              'O. Eosinophil count (individual level, spline w/ 12 knots)\n', 
                              'P. Eosinophil count (individual level, spline w/ 8 knots, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_12_knots_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                             plot_group_mcv_12_knots, plot_group_mcv, plot_ind_mcv_12_knots, plot_ind_mcv, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_mch_12_knots, plot_group_mch, plot_ind_mch_12_knots, plot_ind_mch,
                             NULL, NULL, NULL, NULL, 
                             plot_group_rbc_12_knots, plot_group_rbc, plot_ind_rbc_12_knots, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level, spline w/ 12 knots)\n', 
                              'B. Mean corpuscular volume (group level, spline w/ 8 knots,\nmain text ver.)',
                              'C. Mean corpuscular volume (individual level, spline w/ 12\nknots)',
                              'D. Mean corpuscular volume (individual level, spline w/ 8 knots,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, spline w/ 12\nknots)',
                              'F. Mean corpuscular hemoglobin (group level, spline w/ 8 knots,\nmain text ver.)',
                              'G. Mean corpuscular hemoglobin (individual level, spline w/ 12\nknots)',
                              'H. Mean corpuscular hemoglobin (individual level, spline w/ 8\nknots, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, spline w/ 12 knots)\n', 
                              'J. Red blood cell count (group level, spline w/ 8 knots, main text\nver.)',
                              'K. Red blood cell count (individual level, spline w/ 12 knots)\n', 
                              'L. Red blood cell count (individual level, spline w/ 8 knots, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_12_knots_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                             plot_group_cystatin_c_12_knots, plot_group_cystatin_c, plot_ind_cystatin_c_12_knots, plot_ind_cystatin_c, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_platelet_12_knots, plot_group_platelet, plot_ind_platelet_12_knots, plot_ind_platelet, 
                             NULL, NULL, NULL, NULL, 
                             plot_group_triglycerides_12_knots, plot_group_triglycerides, plot_ind_triglycerides_12_knots, plot_ind_triglycerides,
                             NULL, NULL, NULL, NULL, 
                             plot_group_ldl_12_knots, plot_group_ldl, plot_ind_ldl_12_knots, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, spline w/ 12 knots)\n', 
                              'B. Cystatin C level (group level, spline w/ 8 knots, main text ver.)\n',
                              'C. Cystatin C level (individual level, spline w/ 12 knots)\n', 
                              'D. Cystatin C level (individual level, spline w/ 8 knots, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, spline w/ 12 knots)\n', 
                              'F. Platelet count (group level, spline w/ 8 knots, main text ver.)\n', 
                              'G. Platelet count (individual level, spline w/ 12 knots)\n',
                              'H. Platelet count (individual level, spline w/ 8 knots, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, spline w/ 12 knots)\n', 
                              'J. Triglyceride level (group level, spline w/ 8 knots, main text\nver.)', 
                              'K. Triglyceride level (individual level, spline w/ 12 knots)\n', 
                              'L. Triglyceride level (individual level, spline w/ 8 knots, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, spline w/ 12 knots)\n', 
                              'N. LDL cholesterol level (group level, spline w/ 8 knots, main\ntext ver.)',
                              'O. LDL cholesterol level (individual level, spline w/ 12 knots)\n', 
                              'P. LDL cholesterol level (individual level, spline w/ 8 knots,\nmain text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s46_12_knots_physical.pdf", width = 48, height = 24.82, onefile = T)
grid.arrange(arrangeGrob(plot_12_knots_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s47_12_knots_wbc.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_12_knots_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s48_12_knots_rbc.pdf", width = 48, height = 20.45, onefile = T)
grid.arrange(arrangeGrob(plot_12_knots_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s49_12_knots_other.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_12_knots_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= GWAS with assessment center and genotype array =================================
# Read in group level and individual PGS files for covariates
group_pgs_df_array_center <- read_tsv("data/pgs_pred/group_pgs_df_array_center.tsv") %>%
  filter(threshold == 1)
ind_pgs_df_array_center <- read_tsv("data/pgs_pred/ind_pgs_df_array_center.tsv") %>%
  select(-(16:34))

# Group level plot
plot_group_level_array_center <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
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
      annotate(geom = "segment", x = 103, y = 0.75, xend = 110, yend = 0.75, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 278\nindividuals", x = 130, y = 0.75, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_array_center <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_pred_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_pred_error) %>%
    select(-ref_pred_error) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(pc_dist > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
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
  return(plot)
}

# Make the plots
plot_group_height_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Height", upper)
plot_ind_height_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Height", upper)

plot_group_cystatin_c_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Cystatin_C", upper)
plot_ind_cystatin_c_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Cystatin_C", upper)

plot_group_platelet_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Platelet", upper)
plot_ind_platelet_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Platelet", upper)

plot_group_mcv_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "MCV", upper)
plot_ind_mcv_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "MCV", upper)

plot_group_weight_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Weight", upper)
plot_ind_weight_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Weight", upper)

plot_group_mch_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "MCH", upper)
plot_ind_mch_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "MCH", upper)

plot_group_bmi_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "BMI", upper)
plot_ind_bmi_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "BMI", upper)

plot_group_rbc_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "RBC", upper)
plot_ind_rbc_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "RBC", upper)

plot_group_body_fat_perc_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Body_Fat_Perc", upper)

plot_group_monocyte_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Monocyte", upper)
plot_ind_monocyte_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Monocyte", upper)

plot_group_triglycerides_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Triglycerides", upper)
plot_ind_triglycerides_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Triglycerides", upper)

plot_group_lymphocyte_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Lymphocyte", upper)
plot_ind_lymphocyte_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Lymphocyte", upper)

plot_group_wbc_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "WBC", upper)
plot_ind_wbc_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "WBC", upper)

plot_group_eosinophil_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "Eosinophil", upper)
plot_ind_eosinophil_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "Eosinophil", upper)

plot_group_ldl_array_center <- plot_group_level_array_center(group_pgs_df_array_center, "LDL", upper)
plot_ind_ldl_array_center <- plot_ind_level_array_center(ind_pgs_df_array_center, "LDL", upper)

plot_array_center_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                                 plot_group_height_array_center, plot_group_height, plot_ind_height_array_center, plot_ind_height, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_weight_array_center, plot_group_weight, plot_ind_weight_array_center, plot_ind_weight,
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_bmi_array_center, plot_group_bmi, plot_ind_bmi_array_center, plot_ind_bmi, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_body_fat_perc_array_center, plot_group_body_fat_perc, plot_ind_body_fat_perc_array_center, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, GWAS also controlling for assessment\ncenter and genotype array)',
                              'B. Height (group level, GWAS w/o controlling for assessment\ncenter or genotype array, main text ver.)',
                              'C. Height (individual level, GWAS also controlling for\nassessment center and genotype array)',
                              'D. Height (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, GWAS also controlling for assessment\ncenter and genotype array)',
                              'F. Weight (group level, GWAS w/o controlling for assessment\ncenter or genotype array, main text ver.)',
                              'G. Weight (individual level, GWAS also controlling for\nassessment center and genotype array)',
                              'H. Weight (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, GWAS also controlling for assessment\ncenter and genotype array)',
                              'J. BMI (group level, GWAS w/o controlling for assessment\ncenter or genotype array, main text ver.)',
                              'K. BMI (individual level, GWAS also controlling for\nassessment center and genotype array)',
                              'L. BMI (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'N. Body fat percentage (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'O. Body fat percentage (individual level, GWAS also controlling\nfor assessment center and genotype array)',
                              'P. Body fat percentage (individual level, GWAS w/o controlling\nfor assessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.28, 1, 0.28, 1, 0.28, 1, 0.28, 1),
                   label_fontfamily = "Helvetica")

plot_array_center_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                                 plot_group_monocyte_array_center, plot_group_monocyte, plot_ind_monocyte_array_center, plot_ind_monocyte, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_lymphocyte_array_center, plot_group_lymphocyte, plot_ind_lymphocyte_array_center, plot_ind_lymphocyte, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_wbc_array_center, plot_group_wbc, plot_ind_wbc_array_center, plot_ind_wbc, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_eosinophil_array_center, plot_group_eosinophil, plot_ind_eosinophil_array_center, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'B. Monocyte count (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'C. Monocyte count (individual level, GWAS also controlling for\nassessment center and genotype array)',
                              'D. Monocyte count (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'F. Lymphocyte count (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'G. Lymphocyte count (individual level, GWAS also controlling\nfor assessment center and genotype array)',
                              'H. Lymphocyte count (individual level, GWAS w/o controlling\nfor assessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, GWAS also controlling\nfor assessment center and genotype array)\n',
                              'J. White blood cell count (group level, GWAS w/o controlling\nfor assessment center or genotype array, main text ver.)\n',
                              'K. White blood cell count (individual level, GWAS also\ncontrolling for assessment center and genotype array)\n',
                              'L. White blood cell count (individual level, GWAS w/o\ncontrolling for assessment center or genotype array, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'N. Eosinophil count (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'O. Eosinophil count (individual level, GWAS also controlling\nfor assessment center and genotype array)',
                              'P. Eosinophil count (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.28, 1, 0.28, 1, 0.4, 1, 0.28, 1),
                   label_fontfamily = "Helvetica")

plot_array_center_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                                 plot_group_mcv_array_center, plot_group_mcv, plot_ind_mcv_array_center, plot_ind_mcv, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_mch_array_center, plot_group_mch, plot_ind_mch_array_center, plot_ind_mch,
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_rbc_array_center, plot_group_rbc, plot_ind_rbc_array_center, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level,  GWAS also\ncontrolling for assessment center and genotype array)\n',
                              'B. Mean corpuscular volume (group level, GWAS w/o\ncontrolling for assessment center or genotype array, main text\nver.)',
                              'C. Mean corpuscular volume (individual level, GWAS also\ncontrolling for assessment center and genotype array)\n',
                              'D. Mean corpuscular volume (individual level, GWAS w/o\ncontrolling for assessment center or genotype array, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, GWAS also\ncontrolling for assessment center and genotype array)\n',
                              'F. Mean corpuscular hemoglobin (group level, GWAS w/o\ncontrolling for assessment center or genotype array, main text\nver.)',
                              'G. Mean corpuscular hemoglobin (individual level, GWAS also\ncontrolling for assessment center and genotype array)\n',
                              'H. Mean corpuscular hemoglobin (individual level, GWAS w/o\ncontrolling for assessment center or genotype array, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'J. Red blood cell count (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'K. Red blood cell count (individual level, GWAS also controlling\nfor assessment center and genotype array)',
                              'L. Red blood cell count (individual level, GWAS w/o controlling\nfor assessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.4, 1, 0.4, 1, 0.28, 1),
                   label_fontfamily = "Helvetica")

plot_array_center_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                                 plot_group_cystatin_c_array_center, plot_group_cystatin_c, plot_ind_cystatin_c_array_center, plot_ind_cystatin_c, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_platelet_array_center, plot_group_platelet, plot_ind_platelet_array_center, plot_ind_platelet, 
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_triglycerides_array_center, plot_group_triglycerides, plot_ind_triglycerides_array_center, plot_ind_triglycerides,
                                 NULL, NULL, NULL, NULL, 
                                 plot_group_ldl_array_center, plot_group_ldl, plot_ind_ldl_array_center, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'B. Cystatin C level (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'C. Cystatin C level (individual level, GWAS also controlling for\nassessment center and genotype array)',
                              'D. Cystatin C level (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'F. Platelet count (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'G. Platelet count (individual level, GWAS also controlling\nfor assessment center and genotype array)',
                              'H. Platelet count (individual level, GWAS w/o controlling\nfor assessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'J. Triglyceride level (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'K. Triglyceride level (individual level, GWAS also controlling for\nassessment center and genotype array)',
                              'L. Triglyceride level (individual level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, GWAS also controlling for\nassessment center and genotype array)',
                              'N. LDL cholesterol level (group level, GWAS w/o controlling for\nassessment center or genotype array, main text ver.)',
                              'O. LDL cholesterol level (individual level, GWAS also\ncontrolling for assessment center and genotype array)',
                              'P. LDL cholesterol level (individual level, GWAS w/o controlling\nfor assessment center or genotype array, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.28, 1, 0.28, 1, 0.28, 1, 0.28, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s50_array_center_physical.pdf", width = 48, height = 27.93, onefile = T)
grid.arrange(arrangeGrob(plot_array_center_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s51_array_center_wbc.pdf", width = 48, height = 28.58, onefile = T)
grid.arrange(arrangeGrob(plot_array_center_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s52_array_center_rbc.pdf", width = 48, height = 22.25, onefile = T)
grid.arrange(arrangeGrob(plot_array_center_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s53_array_center_other.pdf", width = 48, height = 27.93, onefile = T)
grid.arrange(arrangeGrob(plot_array_center_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= GWAS with regenie =================================
# Read in group level and individual PGS files for covariates
group_pgs_df_regenie<- read_tsv("data/pgs_pred/group_pgs_df_regenie.tsv") %>%
  filter(threshold == 1)
ind_pgs_df_regenie <- read_tsv("data/pgs_pred/ind_pgs_df_regenie.tsv") %>%
  select(-(16:34))

# Group level plot
plot_group_level_regenie <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
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
      annotate(geom = "segment", x = 96, y = 1.15, xend = 103, yend = 1.15, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 278\nindividuals", x = 123, y = 1.15, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_regenie <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
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
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
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
  return(plot)
}

# Make the plots
plot_group_height_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Height", upper)
plot_ind_height_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Height", upper)

plot_group_cystatin_c_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Cystatin_C", upper)
plot_ind_cystatin_c_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Cystatin_C", upper)

plot_group_platelet_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Platelet", upper)
plot_ind_platelet_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Platelet", upper)

plot_group_mcv_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "MCV", upper)
plot_ind_mcv_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "MCV", upper)

plot_group_weight_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Weight", upper)
plot_ind_weight_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Weight", upper)

plot_group_mch_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "MCH", upper)
plot_ind_mch_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "MCH", upper)

plot_group_bmi_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "BMI", upper)
plot_ind_bmi_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "BMI", upper)

plot_group_rbc_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "RBC", upper)
plot_ind_rbc_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "RBC", upper)

plot_group_body_fat_perc_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Body_Fat_Perc", upper)

plot_group_monocyte_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Monocyte", upper)
plot_ind_monocyte_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Monocyte", upper)

plot_group_triglycerides_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Triglycerides", upper)
plot_ind_triglycerides_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Triglycerides", upper)

plot_group_lymphocyte_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Lymphocyte", upper)
plot_ind_lymphocyte_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Lymphocyte", upper)

plot_group_wbc_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "WBC", upper)
plot_ind_wbc_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "WBC", upper)

plot_group_eosinophil_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "Eosinophil", upper)
plot_ind_eosinophil_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "Eosinophil", upper)

plot_group_ldl_regenie <- plot_group_level_regenie(group_pgs_df_regenie, "LDL", upper)
plot_ind_ldl_regenie <- plot_ind_level_regenie(ind_pgs_df_regenie, "LDL", upper)

plot_regenie_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                            plot_group_height_regenie, plot_group_height, plot_ind_height_regenie, plot_ind_height, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_weight_regenie, plot_group_weight, plot_ind_weight_regenie, plot_ind_weight,
                            NULL, NULL, NULL, NULL, 
                            plot_group_bmi_regenie, plot_group_bmi, plot_ind_bmi_regenie, plot_ind_bmi, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_body_fat_perc_regenie, plot_group_body_fat_perc, plot_ind_body_fat_perc_regenie, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, GWAS w/ regenie)',
                              'B. Height (group level, GWAS w/ PLINK, main text ver.)',
                              'C. Height (individual level, GWAS w/ regenie)',
                              'D. Height (individual level, GWAS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, GWAS w/ regenie)',
                              'F. Weight (group level, GWAS w/ PLINK, main text ver.)',
                              'G. Weight (individual level, GWAS w/ regenie)',
                              'H. Weight (individual level, GWAS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, GWAS w/ regenie)',
                              'J. BMI (group level, GWAS w/ PLINK, main text ver.)',
                              'K. BMI (individual level, GWAS w/ regenie)',
                              'L. BMI (individual level, GWAS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, GWAS w/ regenie)\n',
                              'N. Body fat percentage (group level, GWAS w/ PLINK, main text\nver.)',
                              'O. Body fat percentage (individual level, GWAS w/ regenie)\n',
                              'P. Body fat percentage (individual level, GWAS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_regenie_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                            plot_group_monocyte_regenie, plot_group_monocyte, plot_ind_monocyte_regenie, plot_ind_monocyte, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_lymphocyte_regenie, plot_group_lymphocyte, plot_ind_lymphocyte_regenie, plot_ind_lymphocyte, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_wbc_regenie, plot_group_wbc, plot_ind_wbc_regenie, plot_ind_wbc, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_eosinophil_regenie, plot_group_eosinophil, plot_ind_eosinophil_regenie, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, GWAS w/ regenie)\n',
                              'B. Monocyte count (group level, GWAS w/ PLINK, main text ver.)\n',
                              'C. Monocyte count (individual level, GWAS w/ regenie)\n',
                              'D. Monocyte count (individual level, GWAS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, GWAS w/ regenie)\n',
                              'F. Lymphocyte count (group level, GWAS w/ PLINK, main text\nver.)',
                              'G. Lymphocyte count (individual level, GWAS w/ regenie)\n',
                              'H. Lymphocyte count (individual level, GWAS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, GWAS w/ regenie)\n',
                              'J. White blood cell count (group level, GWAS w/ PLINK, main\ntext ver.)',
                              'K. White blood cell count (individual level, GWAS w/ regenie)\n',
                              'L. White blood cell count (individual level, GWAS w/ PLINK,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, GWAS w/ regenie)\n',
                              'N. Eosinophil count (group level, GWAS w/ PLINK, main text\nver.)',
                              'O. Eosinophil count (individual level, GWAS w/ regenie)\n',
                              'P. Eosinophil count (individual level, GWAS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_regenie_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                            plot_group_mcv_regenie, plot_group_mcv, plot_ind_mcv_regenie, plot_ind_mcv, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_mch_regenie, plot_group_mch, plot_ind_mch_regenie, plot_ind_mch,
                            NULL, NULL, NULL, NULL, 
                            plot_group_rbc_regenie, plot_group_rbc, plot_ind_rbc_regenie, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level, GWAS w/ regenie)\n',
                              'B. Mean corpuscular volume (group level, GWAS w/ PLINK,\nmain text ver.)',
                              'C. Mean corpuscular volume (individual level, GWAS w/\nregenie)',
                              'D. Mean corpuscular volume (individual level, GWAS w/ PLINK,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, GWAS w/\nregenie)',
                              'F. Mean corpuscular hemoglobin (group level, GWAS w/ PLINK,\nmain text ver.)',
                              'G. Mean corpuscular hemoglobin (individual level, GWAS w/\nregenie)',
                              'H. Mean corpuscular hemoglobin (individual level, GWAS w/\nPLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, GWAS w/ regenie)\n',
                              'J. Red blood cell count (group level, GWAS w/ PLINK, main text\nver.)',
                              'K. Red blood cell count (individual level, GWAS w/ regenie)\n',
                              'L. Red blood cell count (individual level, GWAS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_regenie_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                            plot_group_cystatin_c_regenie, plot_group_cystatin_c, plot_ind_cystatin_c_regenie, plot_ind_cystatin_c, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_platelet_regenie, plot_group_platelet, plot_ind_platelet_regenie, plot_ind_platelet, 
                            NULL, NULL, NULL, NULL, 
                            plot_group_triglycerides_regenie, plot_group_triglycerides, plot_ind_triglycerides_regenie, plot_ind_triglycerides,
                            NULL, NULL, NULL, NULL, 
                            plot_group_ldl_regenie, plot_group_ldl, plot_ind_ldl_regenie, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, GWAS w/ regenie)\n',
                              'B. Cystatin C level (group level, GWAS w/ PLINK, main text ver.)\n',
                              'C. Cystatin C level (individual level, GWAS w/ regenie)\n',
                              'D. Cystatin C level (individual level, GWAS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, GWAS w/ regenie)\n',
                              'F. Platelet count (group level, GWAS w/ PLINK, main text ver.)\n',
                              'G. Platelet count (individual level, GWAS w/ regenie)\n',
                              'H. Platelet count (individual level, GWAS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, GWAS w/ regenie)\n',
                              'J. Triglyceride level (group level, GWAS w/ PLINK, main text\nver.)',
                              'K. Triglyceride level (individual level, GWAS w/ regenie)\n',
                              'L. Triglyceride level (individual level, GWAS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, GWAS w/ regenie)\n',
                              'N. LDL cholesterol level (group level, GWAS w/ PLINK, main\ntext ver.)',
                              'O. LDL cholesterol level (individual level, GWAS w/ regenie)\n',
                              'P. LDL cholesterol level (individual level, GWAS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s54_regenie_physical.pdf", width = 48, height = 24.82, onefile = T)
grid.arrange(arrangeGrob(plot_regenie_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s55_regenie_wbc.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_regenie_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s56_regenie_rbc.pdf", width = 48, height = 20.45, onefile = T)
grid.arrange(arrangeGrob(plot_regenie_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s57_regenie_other.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_regenie_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

#================================= PGS with PRS-CS =================================
# Read in group level and individual PGS files for covariates
group_pgs_df_prscs<- read_tsv("data/pgs_pred/group_pgs_df_prscs.tsv")
ind_pgs_df_prscs <- read_tsv("data/pgs_pred/ind_pgs_df_prscs.tsv") %>%
  select(-(16:34))

# Group level plot
plot_group_level_prscs <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(start, stop, by = step))
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
      annotate("text", label = "A bin of 278\nindividuals", x = 135, y = 0.7, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_prscs <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 50 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
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
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
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
  
  if(trait %in% c("MCV", "MCH", "Eosinophil", "Cystatin_C")){
    ymin <- 0.4
    ymax <- 4.2
  } else if(trait %in% c("Lymphocyte", "WBC")){
    ymin <- 0.18
    ymax <- 2.2
  } else if(trait == "Monocyte"){
    ymin <- 0.35
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
  return(plot)
}

# Make the plots
plot_group_height_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Height", upper)
plot_ind_height_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Height", upper)

plot_group_cystatin_c_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Cystatin_C", upper)
plot_ind_cystatin_c_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Cystatin_C", upper)

plot_group_platelet_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Platelet", upper)
plot_ind_platelet_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Platelet", upper)

plot_group_mcv_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "MCV", upper)
plot_ind_mcv_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "MCV", upper)

plot_group_weight_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Weight", upper)
plot_ind_weight_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Weight", upper)

plot_group_mch_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "MCH", upper)
plot_ind_mch_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "MCH", upper)

plot_group_bmi_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "BMI", upper)
plot_ind_bmi_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "BMI", upper)

plot_group_rbc_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "RBC", upper)
plot_ind_rbc_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "RBC", upper)

plot_group_body_fat_perc_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Body_Fat_Perc", upper)

plot_group_monocyte_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Monocyte", upper)
plot_ind_monocyte_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Monocyte", upper)

plot_group_triglycerides_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Triglycerides", upper)
plot_ind_triglycerides_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Triglycerides", upper)

plot_group_lymphocyte_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Lymphocyte", upper)
plot_ind_lymphocyte_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Lymphocyte", upper)

plot_group_wbc_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "WBC", upper)
plot_ind_wbc_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "WBC", upper)

plot_group_eosinophil_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "Eosinophil", upper)
plot_ind_eosinophil_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "Eosinophil", upper)

plot_group_ldl_prscs <- plot_group_level_prscs(group_pgs_df_prscs, "LDL", upper)
plot_ind_ldl_prscs <- plot_ind_level_prscs(ind_pgs_df_prscs, "LDL", upper)

plot_prscs_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_height_prscs, plot_group_height, plot_ind_height_prscs, plot_ind_height, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_weight_prscs, plot_group_weight, plot_ind_weight_prscs, plot_ind_weight,
                          NULL, NULL, NULL, NULL, 
                          plot_group_bmi_prscs, plot_group_bmi, plot_ind_bmi_prscs, plot_ind_bmi, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_body_fat_perc_prscs, plot_group_body_fat_perc, plot_ind_body_fat_perc_prscs, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, PGS w/ PRS-CS)',
                              'B. Height (group level, PGS w/ PLINK, main text ver.)',
                              'C. Height (individual level, PGS w/ PRS-CS)',
                              'D. Height (individual level, PGS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, PGS w/ PRS-CS)',
                              'F. Weight (group level, PGS w/ PLINK, main text ver.)',
                              'G. Weight (individual level, PGS w/ PRS-CS)',
                              'H. Weight (individual level, PGS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, PGS w/ PRS-CS)',
                              'J. BMI (group level, PGS w/ PLINK, main text ver.)',
                              'K. BMI (individual level, PGS w/ PRS-CS)',
                              'L. BMI (individual level, PGS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, PGS w/ PRS-CS)\n',
                              'N. Body fat percentage (group level, PGS w/ PLINK, main text\nver.)',
                              'O. Body fat percentage (individual level, PGS w/ PRS-CS)\n',
                              'P. Body fat percentage (individual level, PGS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_prscs_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_monocyte_prscs, plot_group_monocyte, plot_ind_monocyte_prscs, plot_ind_monocyte, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_lymphocyte_prscs, plot_group_lymphocyte, plot_ind_lymphocyte_prscs, plot_ind_lymphocyte, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_wbc_prscs, plot_group_wbc, plot_ind_wbc_prscs, plot_ind_wbc, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_eosinophil_prscs, plot_group_eosinophil, plot_ind_eosinophil_prscs, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, PGS w/ PRS-CS)\n',
                              'B. Monocyte count (group level, PGS w/ PLINK, main text ver.)\n',
                              'C. Monocyte count (individual level, PGS w/ PRS-CS)\n',
                              'D. Monocyte count (individual level, PGS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, PGS w/ PRS-CS)\n',
                              'F. Lymphocyte count (group level, PGS w/ PLINK, main text\nver.)',
                              'G. Lymphocyte count (individual level, PGS w/ PRS-CS)\n',
                              'H. Lymphocyte count (individual level, PGS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, PGS w/ PRS-CS)\n',
                              'J. White blood cell count (group level, PGS w/ PLINK, main text\nver.)',
                              'K. White blood cell count (individual level, PGS w/ PRS-CS)\n',
                              'L. White blood cell count (individual level, PGS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, PGS w/ PRS-CS)\n',
                              'N. Eosinophil count (group level, PGS w/ PLINK, main text ver.)\n',
                              'O. Eosinophil count (individual level, PGS w/ PRS-CS)\n',
                              'P. Eosinophil count (individual level, PGS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_prscs_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_mcv_prscs, plot_group_mcv, plot_ind_mcv_prscs, plot_ind_mcv, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_mch_prscs, plot_group_mch, plot_ind_mch_prscs, plot_ind_mch,
                          NULL, NULL, NULL, NULL, 
                          plot_group_rbc_prscs, plot_group_rbc, plot_ind_rbc_prscs, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level, PGS w/ PRS-CS)\n',
                              'B. Mean corpuscular volume (group level, PGS w/ PLINK, main\ntext ver.)',
                              'C. Mean corpuscular volume (individual level, PGS w/ PRS-CS)\n',
                              'D. Mean corpuscular volume (individual level, PGS w/ PLINK,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, PGS w/ PRS-CS)\n',
                              'F. Mean corpuscular hemoglobin (group level, PGS w/ PLINK,\nmain text ver.)',
                              'G. Mean corpuscular hemoglobin (individual level, PGS w/\nPRS-CS)',
                              'H. Mean corpuscular hemoglobin (individual level, PGS w/\nPLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, PGS w/ PRS-CS)\n',
                              'J. Red blood cell count (group level, PGS w/ PLINK, main text\nver.)',
                              'K. Red blood cell count (individual level, PGS w/ PRS-CS)\n',
                              'L. Red blood cell count (individual level, PGS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_prscs_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_cystatin_c_prscs, plot_group_cystatin_c, plot_ind_cystatin_c_prscs, plot_ind_cystatin_c, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_platelet_prscs, plot_group_platelet, plot_ind_platelet_prscs, plot_ind_platelet, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_triglycerides_prscs, plot_group_triglycerides, plot_ind_triglycerides_prscs, plot_ind_triglycerides,
                          NULL, NULL, NULL, NULL, 
                          plot_group_ldl_prscs, plot_group_ldl, plot_ind_ldl_prscs, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, PGS w/ PRS-CS)\n',
                              'B. Cystatin C level (group level, PGS w/ PLINK, main text ver.)\n',
                              'C. Cystatin C level (individual level, PGS w/ PRS-CS)\n',
                              'D. Cystatin C level (individual level, PGS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, PGS w/ PRS-CS)',
                              'F. Platelet count (group level, PGS w/ PLINK, main text ver.)',
                              'G. Platelet count (individual level, PGS w/ PRS-CS)',
                              'H. Platelet count (individual level, PGS w/ PLINK, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, PGS w/ PRS-CS)\n',
                              'J. Triglyceride level (group level, PGS w/ PLINK, main text ver.)\n',
                              'K. Triglyceride level (individual level, PGS w/ PRS-CS)\n',
                              'L. Triglyceride level (individual level, PGS w/ PLINK, main text\nver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, PGS w/ PRS-CS)\n',
                              'N. LDL cholesterol level (group level, PGS w/ PLINK, main text\nver.)',
                              'O. LDL cholesterol level (individual level, PGS w/ PRS-CS)\n',
                              'P. LDL cholesterol level (individual level, PGS w/ PLINK, main\ntext ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.1, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s58_prscs_physical.pdf", width = 48, height = 24.82, onefile = T)
grid.arrange(arrangeGrob(plot_prscs_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s59_prscs_wbc.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_prscs_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s60_prscs_rbc.pdf", width = 48, height = 20.45, onefile = T)
grid.arrange(arrangeGrob(plot_prscs_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s61_prscs_other.pdf", width = 48, height = 26.45, onefile = T)
grid.arrange(arrangeGrob(plot_prscs_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= PC distance with 16 PCs =================================
# Read in group level and individual PGS files for covariates
group_pgs_df_16 <- read_tsv("data/pgs_pred/group_pgs_df_16.tsv") %>%
  filter(threshold == 1)
ind_pgs_df_16 <- read_tsv("data/pgs_pred/ind_pgs_df_16.tsv") %>%
  select(-(16:34))

# Upper limit for the GWAS set
pc_dist_gwas_16 <- read_tsv("data/pca/pc_dist_16_gwas_std.tsv")
upper_16 <- unname(quantile(pc_dist_gwas_16$pc_dist, 0.975))

# Get positions of the knots by density
temp_16 <- ind_pgs_df_16 %>% filter(phenotype == "Height" & pc_dist > upper_16)
temp_16 <- temp_16 %>% arrange(pc_dist)
knots_16 <- c(upper_16,
              temp_16$pc_dist[round(nrow(temp_16) / 9)], 
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 2)], 
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 3)], 
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 4)],
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 5)], 
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 6)], 
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 7)], 
              temp_16$pc_dist[round(nrow(temp_16) / 9 * 8)])

# Group level plot
plot_group_level_16 <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots_16), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(upper, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = median_pc, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200), expand = c(0, 0)) +
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
    coord_cartesian(xlim = c(upper, 210), ylim = c(0.0, ymax))
  
  # Add some labels for height
  if(trait == "Height"){
    plot <- plot +
      annotate("text", label = "↑", x = 157, y = 1.694, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 202, y = 1.694, size = 8, hjust = 1, 
               vjust = 0.8, family = "Helvetica") +
      annotate(geom = "segment", x = 116, y = 0.7, xend = 123, yend = 0.7, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 278\nindividuals", x = 143, y = 0.7, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_16 <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_pred_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_pred_error) %>%
    select(-ref_pred_error) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(pc_dist > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots_16), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200), expand = c(0, 0)) +
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
  
  if(trait %in% c("MCV", "MCH", "Eosinophil", "Cystatin_C")){
    ymin <- 0.4
    ymax <- 4.2
  } else if(trait %in% c("Lymphocyte", "WBC")){
    ymin <- 0.18
    ymax <- 2.2
  } else if(trait == "Monocyte"){
    ymin <- 0.35
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
    coord_cartesian(xlim = c(upper, 210), ylim = c(ymax, ymin))
  
  # Add some labels for height
  if(trait == "Height"){
    plot <- plot +
      annotate(geom = "segment", x = 113, y = 0.75, xend = 120, yend = 0.75, 
               color = "black", size = 1) +
      annotate("text", label = "An individual", x = 140, y = 0.75, size = 8, 
               family = "Helvetica") +
      annotate("text", label = "↑", x = 157, y = 0.45, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 202, y = 0.45, size = 8, hjust = 1, vjust = 0.8,
               family = "Helvetica")
  }
  return(plot)
}

# Make the plots
plot_group_height_16 <- plot_group_level_16(group_pgs_df_16, "Height", upper_16)
plot_ind_height_16 <- plot_ind_level_16(ind_pgs_df_16, "Height", upper_16)

plot_group_cystatin_c_16 <- plot_group_level_16(group_pgs_df_16, "Cystatin_C", upper_16)
plot_ind_cystatin_c_16 <- plot_ind_level_16(ind_pgs_df_16, "Cystatin_C", upper_16)

plot_group_platelet_16 <- plot_group_level_16(group_pgs_df_16, "Platelet", upper_16)
plot_ind_platelet_16 <- plot_ind_level_16(ind_pgs_df_16, "Platelet", upper_16)

plot_group_mcv_16 <- plot_group_level_16(group_pgs_df_16, "MCV", upper_16)
plot_ind_mcv_16 <- plot_ind_level_16(ind_pgs_df_16, "MCV", upper_16)

plot_group_weight_16 <- plot_group_level_16(group_pgs_df_16, "Weight", upper_16)
plot_ind_weight_16 <- plot_ind_level_16(ind_pgs_df_16, "Weight", upper_16)

plot_group_mch_16 <- plot_group_level_16(group_pgs_df_16, "MCH", upper_16)
plot_ind_mch_16 <- plot_ind_level_16(ind_pgs_df_16, "MCH", upper_16)

plot_group_bmi_16 <- plot_group_level_16(group_pgs_df_16, "BMI", upper_16)
plot_ind_bmi_16 <- plot_ind_level_16(ind_pgs_df_16, "BMI", upper_16)

plot_group_rbc_16 <- plot_group_level_16(group_pgs_df_16, "RBC", upper_16)
plot_ind_rbc_16 <- plot_ind_level_16(ind_pgs_df_16, "RBC", upper_16)

plot_group_body_fat_perc_16 <- plot_group_level_16(group_pgs_df_16, "Body_Fat_Perc", upper_16)
plot_ind_body_fat_perc_16 <- plot_ind_level_16(ind_pgs_df_16, "Body_Fat_Perc", upper_16)

plot_group_monocyte_16 <- plot_group_level_16(group_pgs_df_16, "Monocyte", upper_16)
plot_ind_monocyte_16 <- plot_ind_level_16(ind_pgs_df_16, "Monocyte", upper_16)

plot_group_triglycerides_16 <- plot_group_level_16(group_pgs_df_16, "Triglycerides", upper_16)
plot_ind_triglycerides_16 <- plot_ind_level_16(ind_pgs_df_16, "Triglycerides", upper_16)

plot_group_lymphocyte_16 <- plot_group_level_16(group_pgs_df_16, "Lymphocyte", upper_16)
plot_ind_lymphocyte_16 <- plot_ind_level_16(ind_pgs_df_16, "Lymphocyte", upper_16)

plot_group_wbc_16 <- plot_group_level_16(group_pgs_df_16, "WBC", upper_16)
plot_ind_wbc_16 <- plot_ind_level_16(ind_pgs_df_16, "WBC", upper_16)

plot_group_eosinophil_16 <- plot_group_level_16(group_pgs_df_16, "Eosinophil", upper_16)
plot_ind_eosinophil_16 <- plot_ind_level_16(ind_pgs_df_16, "Eosinophil", upper_16)

plot_group_ldl_16 <- plot_group_level_16(group_pgs_df_16, "LDL", upper_16)
plot_ind_ldl_16 <- plot_ind_level_16(ind_pgs_df_16, "LDL", upper_16)

plot_16_pc_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_height_16, plot_group_height, plot_ind_height_16, plot_ind_height, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_weight_16, plot_group_weight, plot_ind_weight_16, plot_ind_weight,
                          NULL, NULL, NULL, NULL, 
                          plot_group_bmi_16, plot_group_bmi, plot_ind_bmi_16, plot_ind_bmi, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_body_fat_perc_16, plot_group_body_fat_perc, plot_ind_body_fat_perc_16, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, genetic distance w/ 16 PCs)\n',
                              'B. Height (group level, genetic distance w/ 40 PCs, main text\nver.)',
                              'C. Height (individual level, genetic distance w/ 16 PCs)\n',
                              'D. Height (individual level, genetic distance w/ 40 PCs, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, genetic distance w/ 16 PCs)\n',
                              'F. Weight (group level, genetic distance w/ 40 PCs, main text\nver.)',
                              'G. Weight (individual level, genetic distance w/ 16 PCs)\n',
                              'H. Weight (individual level, genetic distance w/ 40 PCs, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, genetic distance w/ 16 PCs)\n', 
                              'J. BMI (group level, genetic distance w/ 40 PCs, main text\nver.)',
                              'K. BMI (individual level, genetic distance w/ 16 PCs)\n',
                              'L. BMI (individual level, genetic distance w/ 40 PCs, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, genetic distance w/ 16\nPCs)', 
                              'N. Body fat percentage (group level, genetic distance w/ 40\nPCs, main text ver.)', 
                              'O. Body fat percentage (individual level, genetic distance w/ 16\nPCs)',
                              'P. Body fat percentage (individual level, genetic distance w/ 40\nPCs, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_16_pc_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_monocyte_16, plot_group_monocyte, plot_ind_monocyte_16, plot_ind_monocyte, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_lymphocyte_16, plot_group_lymphocyte, plot_ind_lymphocyte_16, plot_ind_lymphocyte, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_wbc_16, plot_group_wbc, plot_ind_wbc_16, plot_ind_wbc, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_eosinophil_16, plot_group_eosinophil, plot_ind_eosinophil_16, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, genetic distance w/ 16 PCs)\n', 
                              'B. Monocyte count (group level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              'C. Monocyte count (individual level, genetic distance w/ 16 PCs)\n',
                              'D. Monocyte count (individual level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, genetic distance w/ 16 PCs)\n', 
                              'F. Lymphocyte count (group level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              'G. Lymphocyte count (individual level, genetic distance w/ 16\nPCs)', 
                              'H. Lymphocyte count (individual level, genetic distance w/ 40\nPCs, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, genetic distance w/ 16\nPCs)', 
                              'J. White blood cell count (group level, genetic distance w/ 40\nPCs, main text ver.)',
                              'K. White blood cell count (individual level, genetic distance w/\n16 PCs)', 
                              'L. White blood cell count (individual level, genetic distance w/\n40 PCs, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, genetic distance w/ 16 PCs)\n', 
                              'N. Eosinophil count (group level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              'O. Eosinophil count (individual level, genetic distance w/ 16\nPCs)', 
                              'P. Eosinophil count (individual level, genetic distance w/ 40\nPCs, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_16_pc_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_mcv_16, plot_group_mcv, plot_ind_mcv_16, plot_ind_mcv, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_mch_16, plot_group_mch, plot_ind_mch_16, plot_ind_mch,
                          NULL, NULL, NULL, NULL, 
                          plot_group_rbc_16, plot_group_rbc, plot_ind_rbc_16, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level, genetic distance w/\n16 PCs)', 
                              'B. Mean corpuscular volume (group level, genetic distance w/\n40 PCs, main text ver.)',
                              'C. Mean corpuscular volume (individual level, genetic distance\nw/ 16 PCs)',
                              'D. Mean corpuscular volume (individual level, genetic distance\nw/ 40 PCs, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, genetic distance\nw/ 16 PCs)',
                              'F. Mean corpuscular hemoglobin (group level, genetic distance\nw/ 40 PCs, main text ver.)',
                              'G. Mean corpuscular hemoglobin (individual level, genetic\ndistance w/ 16 PCs)',
                              'H. Mean corpuscular hemoglobin (individual level, genetic\ndistance w/ 40 PCs, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, genetic distance w/ 16 PCs)\n', 
                              'J. Red blood cell count (group level, genetic distance w/ 40\nPCs, main text ver.)',
                              'K. Red blood cell count (individual level, genetic distance w/ 16\nPCs)', 
                              'L. Red blood cell count (individual level, genetic distance w/ 40\nPCs, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_16_pc_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_cystatin_c_16, plot_group_cystatin_c, plot_ind_cystatin_c_16, plot_ind_cystatin_c, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_platelet_16, plot_group_platelet, plot_ind_platelet_16, plot_ind_platelet, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_triglycerides_16, plot_group_triglycerides, plot_ind_triglycerides_16, plot_ind_triglycerides,
                          NULL, NULL, NULL, NULL, 
                          plot_group_ldl_16, plot_group_ldl, plot_ind_ldl_16, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, genetic distance w/ 16 PCs)\n', 
                              'B. Cystatin C level (group level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              'C. Cystatin C level (individual level, genetic distance w/ 16 PCs)\n', 
                              'D. Cystatin C level (individual level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, genetic distance w/ 16 PCs)\n', 
                              'F. Platelet count (group level, genetic distance w/ 40 PCs, main\ntext ver.)', 
                              'G. Platelet count (individual level, genetic distance w/ 16 PCs)\n',
                              'H. Platelet count (individual level, genetic distance w/ 40 PCs,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, genetic distance w/ 16 PCs)\n', 
                              'J. Triglyceride level (group level, genetic distance w/ 40 PCs,\nmain text ver.)', 
                              'K. Triglyceride level (individual level, genetic distance w/ 16\nPCs)', 
                              'L. Triglyceride level (individual level, genetic distance w/ 40\nPCs, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, genetic distance w/ 16\nPCs)', 
                              'N. LDL cholesterol level (group level, genetic distance w/ 40\nPCs, main text ver.)',
                              'O. LDL cholesterol level (individual level, genetic distance w/\n16 PCs)', 
                              'P. LDL cholesterol level (individual level, genetic distance w/ 40\nPCs, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s62_16_pcs_physical.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_16_pc_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s63_16_pcs_wbc.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_16_pc_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s64_16_pcs_rbc.pdf", width = 48, height = 20.45, onefile = T)
grid.arrange(arrangeGrob(plot_16_pc_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s65_16_pcs_other.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_16_pc_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= Stratifying by Townsend =================================
# Read in group and individual level PGS files for covariates
group_pgs_df_townsend_2 <- read_tsv("data/pgs_pred/group_non_pgs_df_townsend_2_strata.tsv") %>%
  filter(threshold == 1)
ind_pgs_df_townsend_2 <- read_tsv("data/pgs_pred/ind_pgs_df_townsend_2.tsv") %>%
  select(-(16:34))

# Get positions of the knots by density
temp_townsend_2_higher <- ind_pgs_df_townsend_2 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 2 & pc_dist > upper)
temp_townsend_2_higher = temp_townsend_2_higher %>% arrange(pc_dist)
knots_2_higher <- c(upper,
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9)], 
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 2)], 
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 3)], 
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 4)],
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 5)], 
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 6)], 
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 7)], 
              temp_townsend_2_higher$pc_dist[round(nrow(temp_townsend_2_higher) / 9 * 8)])

temp_townsend_2_lower <- ind_pgs_df_townsend_2 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 1 & pc_dist > upper)
temp_townsend_2_lower = temp_townsend_2_lower %>% arrange(pc_dist)
knots_2_lower <- c(upper,
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9)], 
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 2)], 
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 3)], 
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 4)],
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 5)], 
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 6)], 
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 7)], 
              temp_townsend_2_lower$pc_dist[round(nrow(temp_townsend_2_lower) / 9 * 8)])

# Group level plot
plot_group_level_townsend_2 <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  plot_df$townsend_rank <- ifelse(plot_df$townsend_rank == 2, "Top 50%", "Bottom 50%")
  plot_df$townsend_rank <- factor(plot_df$townsend_rank, levels = c("Top 50%", "Bottom 50%"))
  
  plot_df <- plot_df %>% filter(phenotype == trait & median_values > upper)
  plot_df_lower <- plot_df %>% filter(townsend_rank == "Bottom 50%")
  plot_df_higher <- plot_df %>% filter(townsend_rank == "Top 50%")

  start <- upper
  stop_lower <- max(plot_df_lower$median_pc)
  step_lower <- (stop_lower - start) / 400
  stop_higher <- max(plot_df_higher$median_pc)
  step_higher <- (stop_higher - start) / 400

  # Fit a spline
  lm_lower <- lm(relative_performance ~ 
                  splines::bs(median_values, knots = knots_2_lower), data = plot_df_lower)
  plot_df_2_lower <- cbind.data.frame(median_values = seq(upper, stop_lower, by = step_lower))
  plot_df_2_lower$fitted_val <- unname(predict(lm_lower, newdata = plot_df_2_lower))
  
  lm_higher <- lm(relative_performance ~ 
                  splines::bs(median_values, knots = knots_2_higher), data = plot_df_higher)
  plot_df_2_higher <- cbind.data.frame(median_values = seq(upper, stop_higher, by = step_higher))
  plot_df_2_higher$fitted_val <- unname(predict(lm_higher, newdata = plot_df_2_higher))

  my_cols <- c("#BD0026", "#FED25D")
  plot <- ggplot(plot_df, aes(x = median_pc, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size = 5,alpha = 0.4, fill = townsend_rank, color = townsend_rank, shape = 23)+
    scale_fill_manual(name = "Townsend Deprivation\nIndex", values=my_cols) +
    scale_color_manual(name = "Townsend Deprivation\nIndex", values=my_cols) +
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2_lower,
                  aes(x = median_values, y = fitted_val), size=2.5, color = "#E5BD54") +
    geom_line(data = plot_df_2_higher,
                  aes(x = median_values, y = fitted_val), size=2.5, color = "#AA0022")
  
  # Set the upper limit of the y-axis
  if(trait == "Height"){
    ymax <- 1.8
  }  else if(trait %in% c("Body_Fat_Perc", "WBC")){
    ymax <- 3.2
  } else if(trait %in% c("Weight", "BMI")){
    ymax <- 3.4
  } else if(trait %in% c("Eosinophil")){
    ymax <- 3.8
  } else if(trait %in% c("MCH", "MCV")){
    ymax <- 1.8
  } else if(trait %in% c("Triglycerides", "Cystatin_C", "Monocyte")){
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
      annotate("text", label = "↑", x = 62, y = 1.694, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 62, y = 1.694, size = 8, hjust = 1, 
               vjust = 0.8, family = "Helvetica") +
      annotate(geom = "segment", x = 115, y = 0.8, xend = 122, yend = 0.8, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 278\nindividuals", x = 142, y = 0.8, size = 8, 
               family = "Helvetica") +
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = 1.7, size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 50%", 
               x = 195, y = 1.6, size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Top 50%", 
               x = 195, y = 1.5, size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
  } else if(trait == "Monocyte"){
    plot = plot + 
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = 2.078, size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 50%", 
               x = 195, y = 1.956, size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Top 50%", 
               x = 195, y = 1.833, size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
  } else if(trait == "MCV"){
    plot = plot + 
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = 1.7, size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 50%", 
               x = 195, y = 1.6, size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Top 50%", 
               x = 195, y = 1.5, size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
  } else if(trait == "Cystatin_C"){
    plot = plot + 
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = 2.644, size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 50%", 
               x = 195, y = 2.489, size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Top 50%", 
               x = 195, y = 2.333, size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
  } else{
    plot = plot + theme(legend.position="none")
  }
  return(plot)
}

# Read in individual level PGS file for covariates
ind_pgs_df_townsend_5 <- read_tsv("data/pgs_pred/ind_pgs_df_townsend_5.tsv") %>%
  select(-(16:34))
           
# Get positions of the knots by density
temp_townsend_5_highest <- ind_pgs_df_townsend_5 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 5 & pc_dist > upper)
temp_townsend_5_highest = temp_townsend_5_highest %>% arrange(pc_dist)
knots_5_highest <- c(upper,
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9)], 
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 2)], 
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 3)], 
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 4)],
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 5)], 
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 6)], 
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 7)], 
              temp_townsend_5_highest$pc_dist[round(nrow(temp_townsend_5_highest) / 9 * 8)])

temp_townsend_5_higher <- ind_pgs_df_townsend_5 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 4 & pc_dist > upper)
temp_townsend_5_higher = temp_townsend_5_higher %>% arrange(pc_dist)
knots_5_higher <- c(upper,
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9)], 
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 2)], 
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 3)], 
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 4)],
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 5)], 
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 6)], 
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 7)], 
              temp_townsend_5_higher$pc_dist[round(nrow(temp_townsend_5_higher) / 9 * 8)])

temp_townsend_5_middle <- ind_pgs_df_townsend_5 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 3 & pc_dist > upper)
temp_townsend_5_middle = temp_townsend_5_middle %>% arrange(pc_dist)
knots_5_middle <- c(upper,
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9)], 
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 2)], 
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 3)], 
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 4)],
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 5)], 
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 6)], 
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 7)], 
              temp_townsend_5_middle$pc_dist[round(nrow(temp_townsend_5_middle) / 9 * 8)])

temp_townsend_5_lower <- ind_pgs_df_townsend_5 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 2 & pc_dist > upper)
temp_townsend_5_lower = temp_townsend_5_lower %>% arrange(pc_dist)
knots_5_lower <- c(upper,
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9)], 
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 2)], 
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 3)], 
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 4)],
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 5)], 
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 6)], 
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 7)], 
              temp_townsend_5_lower$pc_dist[round(nrow(temp_townsend_5_lower) / 9 * 8)])

temp_townsend_5_lowest <- ind_pgs_df_townsend_5 %>% filter(phenotype == "Height") %>%
  filter(townsend_rank == 1 & pc_dist > upper)
temp_townsend_5_lowest = temp_townsend_5_lowest %>% arrange(pc_dist)
knots_5_lowest <- c(upper,
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9)], 
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 2)], 
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 3)], 
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 4)],
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 5)], 
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 6)], 
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 7)], 
              temp_townsend_5_lowest$pc_dist[round(nrow(temp_townsend_5_lowest) / 9 * 8)])

# Individual level plot
plot_ind_level_townsend_5 <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 25) %>%
    dplyr::summarise(ref_pred_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_pred_error) %>%
    select(-ref_pred_error) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(pc_dist > upper)

  plot_df$townsend_rank <- ifelse(plot_df$townsend_rank == 5, "Top 20%", 
                                 ifelse(plot_df$townsend_rank == 4, "Top 20-40%",
                                        ifelse(plot_df$townsend_rank == 3, "Middle 20%",
                                               ifelse(plot_df$townsend_rank == 2, "Bottom 20-40%",
                                                      "Bottom 20%"))))
  plot_df$townsend_rank <- factor(plot_df$townsend_rank, levels = c("Top 20%", "Top 20-40%",
                                                                   "Middle 20%", "Bottom 20-40%",
                                                                   "Bottom 20%"))
  
  plot_df = plot_df %>% filter(phenotype == trait & pc_dist > upper)
  plot_df_lowest = plot_df %>% filter(townsend_rank == "Bottom 20%")
  plot_df_lower = plot_df %>% filter(townsend_rank == "Bottom 20-40%")
  plot_df_middle = plot_df %>% filter(townsend_rank == "Middle 20%")
  plot_df_higher = plot_df %>% filter(townsend_rank == "Top 20-40%")
  plot_df_highest = plot_df %>% filter(townsend_rank == "Top 20%")
  
  
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots_16), data = plot_df)

  start <- upper
  stop_highest <- max(plot_df_highest$pc_dist)
  step_highest <- (stop_highest - start) / 400
  stop_higher <- max(plot_df_higher$pc_dist)
  step_higher <- (stop_higher - start) / 400
  stop_middle <- max(plot_df_middle$pc_dist)
  step_middle <- (stop_middle - start) / 400
  stop_lower <- max(plot_df_lower$pc_dist)
  step_lower <- (stop_lower - start) / 400
  stop_lowest <- max(plot_df_lowest$pc_dist)
  step_lowest <- (stop_lowest - start) / 400

  # Fit a spline
  lm_lowest <- lm(relative_performance ~ 
                  splines::bs(pc_dist, knots = knots_5_lowest), data = plot_df_lowest)
  plot_df_2_lowest <- cbind.data.frame(pc_dist = seq(start, stop_lowest, by = step_lowest))
  plot_df_2_lowest$fitted_val <- unname(predict(lm_lowest, newdata = plot_df_2_lowest))
  
  lm_lower <- lm(relative_performance ~ 
                  splines::bs(pc_dist, knots = knots_5_lower), data = plot_df_lower)
  plot_df_2_lower <- cbind.data.frame(pc_dist = seq(start, stop_lower, by = step_lower))
  plot_df_2_lower$fitted_val <- unname(predict(lm_lower, newdata = plot_df_2_lower))
  
  lm_middle <- lm(relative_performance ~ 
                  splines::bs(pc_dist, knots = knots_5_middle), data = plot_df_middle)
  plot_df_2_middle <- cbind.data.frame(pc_dist = seq(start, stop_middle, by = step_middle))
  plot_df_2_middle$fitted_val <- unname(predict(lm_middle, newdata = plot_df_2_middle))
  
  lm_higher <- lm(relative_performance ~ 
                   splines::bs(pc_dist, knots = knots_5_higher), data = plot_df_higher)
  plot_df_2_higher <- cbind.data.frame(pc_dist = seq(start, stop_higher, by = step_higher))
  plot_df_2_higher$fitted_val <- unname(predict(lm_higher, newdata = plot_df_2_higher))
  
  lm_highest <- lm(relative_performance ~ 
                   splines::bs(pc_dist, knots = knots_5_highest), data = plot_df_highest)
  plot_df_2_highest <- cbind.data.frame(pc_dist = seq(start, stop_highest, by = step_highest))
  plot_df_2_highest$fitted_val <- unname(predict(lm_highest, newdata = plot_df_2_highest))

 # Replace with nonzero values for plotting
  plot_df_2_lowest$fitted_val[plot_df_2_lowest$fitted_val <= 0] <- 1e-6
  plot_df_2_lower$fitted_val[plot_df_2_lower$fitted_val <= 0] <- 1e-6
  plot_df_2_middle$fitted_val[plot_df_2_middle$fitted_val <= 0] <- 1e-6
  plot_df_2_higher$fitted_val[plot_df_2_higher$fitted_val <= 0] <- 1e-6
  plot_df_2_highest$fitted_val[plot_df_2_highest$fitted_val <= 0] <- 1e-6

  my_cols = c("#BD0026", "#F03B20", "#FD8D3C", "#FEB24C", "#FED25D")
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.05, shape = 16, aes(fill = townsend_rank, color = townsend_rank)) +
    scale_fill_manual(name = "Townsend Deprivation\nIndex", values=my_cols) +
    scale_color_manual(name = "Townsend Deprivation\nIndex", values=my_cols) +
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2_lowest,
                  aes(x = pc_dist, y = fitted_val), size=2.5, color = "#E5BD54") +
    geom_line(data = plot_df_2_lower,
                  aes(x = pc_dist, y = fitted_val), size=2.5, color = "#E5A044") +
    geom_line(data = plot_df_2_middle,
                  aes(x = pc_dist, y = fitted_val), size=2.5, color = "#E47F36") +
    geom_line(data = plot_df_2_higher,
                  aes(x = pc_dist, y = fitted_val), size=2.5, color = "#D6351D") +
    geom_line(data = plot_df_2_highest,
                  aes(x = pc_dist, y = fitted_val), size=2.5, color = "#AA0022")
  
  plot <- plot + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    theme(axis.title = element_text(size=24, family = "Helvetica"),
          axis.title.x = element_blank(),
          axis.text = element_text(size=20, family = "Helvetica")) +
    ylab("Squared prediction error")
  
  if(trait %in% c("MCV", "MCH", "RBC")){
    ymin <- 0.4
    ymax <- 4.2
  } else if(trait %in% c("WBC", "Triglycerides")){
    ymin <- 0.18
    ymax <- 2.2
  } else if(trait %in% c("Lymphocyte")){
    ymin <- 2.4
    ymax <- 1e-6
  } else if(trait == "Cystatin_C"){
    ymin <- 4.2
    ymax <- 1e-6
  } else if(trait == "Monocyte"){
    ymin <- 2.2
    ymax <- 1e-6
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
    pos <- 2 ^ (log2(0.45) + (0:(6 - 1)) * (log2(0.5) - log2(0.45)))
    
    plot <- plot +
      annotate(geom = "segment", x = 108, y = 0.75, xend = 115, yend = 0.75, 
               color = "black", size = 1) +
      annotate("text", label = "An individual", x = 135, y = 0.75, size = 8, 
               family = "Helvetica") +
      annotate("text", label = "↑", x = 20, y = 0.45, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 62, y = 0.45, size = 8, hjust = 1, vjust = 0.8,
               family = "Helvetica") +
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = pos[1], size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 20%", 
               x = 195, y = pos[2], size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Bottom 20-40%", 
               x = 195, y = pos[3], size = 8, hjust = 1, family = "Helvetica", color = "#E5A044") +
      annotate("text", label = "Middle 20%", 
               x = 195, y = pos[4], size = 8, hjust = 1, family = "Helvetica", color = "#E47F36") +
      annotate("text", label = "Top 20-40%", 
               x = 195, y = pos[5], size = 8, hjust = 1, family = "Helvetica", color = "#D6351D") +
      annotate("text", label = "Top 20%", 
               x = 195, y = pos[6], size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
  }else if(trait == "Monocyte"){
    
    pos <- c(2.742871e-06, 6.7638421e-06, 1.6679443e-05, 
             4.1131033e-05, 0.00010142796, 0.00025011847)
    
    plot <- plot + 
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = pos[1], size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 20%", 
               x = 195, y = pos[2], size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Bottom 20-40%", 
               x = 195, y = pos[3], size = 8, hjust = 1, family = "Helvetica", color = "#E5A044") +
      annotate("text", label = "Middle 20%", 
               x = 195, y = pos[4], size = 8, hjust = 1, family = "Helvetica", color = "#E47F36") +
      annotate("text", label = "Top 20-40%", 
               x = 195, y = pos[5], size = 8, hjust = 1, family = "Helvetica", color = "#D6351D") +
      annotate("text", label = "Top 20%", 
               x = 195, y = pos[6], size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
  } else if(trait == "MCV"){
    
    pos <- c(0.45, 0.5071, 0.5720, 0.6460, 0.7300, 0.8269)
    
    plot <- plot + 
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = pos[1], size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 20%", 
               x = 195, y = pos[2], size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Bottom 20-40%", 
               x = 195, y = pos[3], size = 8, hjust = 1, family = "Helvetica", color = "#E5A044") +
      annotate("text", label = "Middle 20%", 
               x = 195, y = pos[4], size = 8, hjust = 1, family = "Helvetica", color = "#E47F36") +
      annotate("text", label = "Top 20-40%", 
               x = 195, y = pos[5], size = 8, hjust = 1, family = "Helvetica", color = "#D6351D") +
      annotate("text", label = "Top 20%", 
               x = 195, y = pos[6], size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
    
  } else if(trait == "Cystatin_C"){
    
    pos <- c(2.8681907e-06, 7.3612633e-06, 1.8892815e-05, 
             4.8488752e-05, 0.00012444726, 0.00031939614)
    
    plot <- plot + 
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = pos[1], size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 20%", 
               x = 195, y = pos[2], size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Bottom 20-40%", 
               x = 195, y = pos[3], size = 8, hjust = 1, family = "Helvetica", color = "#E5A044") +
      annotate("text", label = "Middle 20%", 
               x = 195, y = pos[4], size = 8, hjust = 1, family = "Helvetica", color = "#E47F36") +
      annotate("text", label = "Top 20-40%", 
               x = 195, y = pos[5], size = 8, hjust = 1, family = "Helvetica", color = "#D6351D") +
      annotate("text", label = "Top 20%", 
               x = 195, y = pos[6], size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
    
  } else if(trait == "BMI"){
    
    pos <- 2 ^ (log2(0.45) + (0:(6 - 1)) * (log2(0.5) - log2(0.45)))
    
    plot <- plot +
      annotate(geom = "segment", x = 108, y = 1.5, xend = 115, yend = 1.5, 
               color = "black", size = 1) +
      annotate("text", label = "An individual", x = 135, y = 1.5, size = 8, 
               family = "Helvetica") +
      annotate("text", label = "↑", x = 20, y = 0.45, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 62, y = 0.45, size = 8, hjust = 1, vjust = 0.8,
               family = "Helvetica") +
      annotate("text", label = "Townsend Deprivation Index", 
               x = 195, y = pos[1], size = 8, hjust = 1, family = "Helvetica") +
      annotate("text", label = "Bottom 20%", 
               x = 195, y = pos[2], size = 8, hjust = 1, family = "Helvetica", color = "#E5BD54") +
      annotate("text", label = "Bottom 20-40%", 
               x = 195, y = pos[3], size = 8, hjust = 1, family = "Helvetica", color = "#E5A044") +
      annotate("text", label = "Middle 20%", 
               x = 195, y = pos[4], size = 8, hjust = 1, family = "Helvetica", color = "#E47F36") +
      annotate("text", label = "Top 20-40%", 
               x = 195, y = pos[5], size = 8, hjust = 1, family = "Helvetica", color = "#D6351D") +
      annotate("text", label = "Top 20%", 
               x = 195, y = pos[6], size = 8, hjust = 1, family = "Helvetica", color = "#AA0022") +
      theme(legend.position="none")
    
  } else {
    plot <- plot + theme(legend.position="none")
  }
  return(plot)
}

# Make the plots
plot_group_height_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Height", upper)
plot_ind_height_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Height")

plot_group_cystatin_c_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Cystatin_C", upper)
plot_ind_cystatin_c_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Cystatin_C")

plot_group_platelet_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Platelet", upper)
plot_ind_platelet_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Platelet", upper)

plot_group_mcv_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "MCV", upper)
plot_ind_mcv_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "MCV", upper)

plot_group_weight_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Weight", upper)
plot_ind_weight_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Weight", upper)

plot_group_mch_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "MCH", upper)
plot_ind_mch_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "MCH", upper)

plot_group_bmi_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "BMI", upper)
plot_ind_bmi_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "BMI", upper)

plot_group_rbc_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "RBC", upper)
plot_ind_rbc_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "RBC", upper)

plot_group_body_fat_perc_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Body_Fat_Perc", upper)
plot_ind_body_fat_perc_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Body_Fat_Perc", upper)

plot_group_monocyte_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Monocyte", upper)
plot_ind_monocyte_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Monocyte", upper)

plot_group_triglycerides_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Triglycerides")
plot_ind_triglycerides_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Triglycerides", upper)

plot_group_lymphocyte_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Lymphocyte")
plot_ind_lymphocyte_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Lymphocyte", upper)

plot_group_wbc_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "WBC")
plot_ind_wbc_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "WBC", upper)

plot_group_eosinophil_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "Eosinophil")
plot_ind_eosinophil_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "Eosinophil", upper)

plot_group_ldl_townsend_2 <- plot_group_level_townsend_2(group_pgs_df_townsend_2, "LDL")
plot_ind_ldl_townsend_5 <- plot_ind_level_townsend_5(ind_pgs_df_townsend_5, "LDL", upper)

plot_townsend_1 <- plot_grid(NULL, NULL, 
                            plot_group_height_townsend_2, plot_ind_height_townsend_5, 
                            NULL, NULL,
                            plot_group_weight_townsend_2, plot_ind_weight_townsend_5,
                            NULL, NULL,
                            plot_group_bmi_townsend_2, NULL, 
                            NULL, NULL, 
                            plot_group_body_fat_perc_townsend_2, plot_ind_body_fat_perc_townsend_5,
                   labels = c('A. Height (group level, stratified by SES)',
                              'B. Height (individual level, stratified by SES)',
                              '',
                              '',
                              'C. Weight (group level, stratified by SES)',
                              'D. Weight (individual level, stratified by SES)',
                              '',
                              '',
                              'E. BMI (group level, stratified by SES)',
                              '',
                              '',
                              '',
                              'F. Body fat percentage (group level, stratified by SES)',
                              'G. Body fat percentage (individual level, stratified by SES)',
                              '',
                              ''), 
                   ncol = 2, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

plot_townsend_2 <- plot_grid(NULL, NULL, 
                            plot_group_monocyte_townsend_2, plot_ind_monocyte_townsend_5, 
                            NULL, NULL,
                            plot_group_lymphocyte_townsend_2, plot_ind_lymphocyte_townsend_5,
                            NULL, NULL,
                            plot_group_wbc_townsend_2, plot_ind_wbc_townsend_5, 
                            NULL, NULL, 
                            plot_group_eosinophil_townsend_2, plot_ind_eosinophil_townsend_5,
                   labels = c('A. Monocyte count (group level, stratified by SES)',
                              'B. Monocyte count (individual level, stratified by SES)',
                              '',
                              '',
                              'C. Lymphocyte count (group level, stratified by SES)',
                              'D. Lymphocyte count (individual level, stratified by SES)',
                              '',
                              '',
                              'E. White blood cell count (group level, stratified by SES)',
                              'F. White blood cell count (individual level, stratified by SES)',
                              '',
                              '',
                              'G. Eosinophil count (group level, stratified by SES)',
                              'H. Eosinophil count (individual level, stratified by SES)',
                              '',
                              ''), 
                   ncol = 2, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

plot_townsend_3 <- plot_grid(NULL, NULL, 
                            plot_group_mcv_townsend_2, plot_ind_mcv_townsend_5, 
                            NULL, NULL,
                            plot_group_mch_townsend_2, plot_ind_mch_townsend_5,
                            NULL, NULL,
                            plot_group_rbc_townsend_2, plot_ind_rbc_townsend_5, 
                   labels = c('A. Mean corpuscular volume (group level, stratified by SES)\n',
                              'B. Mean corpuscular volume (individual level, stratified by\nSES)',
                              '',
                              '',
                              'C. Mean corpuscular hemoglobin (group level, stratified by\nSES)',
                              'D. Mean corpuscular hemoglobin (individual level, stratified\nby SES)',
                              '',
                              '',
                              'E. Red blood cell count (group level, stratified by SES)',
                              'F. Red blood cell count (individual level, stratified by SES)',
                              '',
                              ''), 
                   ncol = 2, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

plot_townsend_4 <- plot_grid(NULL, NULL, 
                            plot_group_cystatin_c_townsend_2, plot_ind_cystatin_c_townsend_5, 
                            NULL, NULL,
                            plot_group_platelet_townsend_2, plot_ind_platelet_townsend_5,
                            NULL, NULL,
                            plot_group_triglycerides_townsend_2, plot_ind_triglycerides_townsend_5, 
                            NULL, NULL, 
                            plot_group_ldl_townsend_2, plot_ind_ldl_townsend_5,
                   labels = c('A. Cystatin C level (group level, stratified by SES)',
                              'B. Cystatin C level (individual level, stratified by SES)',
                              '',
                              '',
                              'C. Platelet count (group level, stratified by SES)',
                              'D. Platelet count (individual level, stratified by SES)',
                              '',
                              '',
                              'E. Triglyceride level (group level, stratified by SES)',
                              'F. Triglyceride level (individual level, stratified by SES)',
                              '',
                              '',
                              'G. LDL cholesterol level (group level, stratified by SES)',
                              'H. LDL cholesterol level (individual level, stratified by SES)',
                              '',
                              ''), 
                   ncol = 2, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s67_townsend_port_physical.pdf", width = 24, height = 24, onefile = T)
grid.arrange(arrangeGrob(plot_townsend_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s68_townsend_port_wbc.pdf", width = 24, height = 24, onefile = T)
grid.arrange(arrangeGrob(plot_townsend_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s69_townsend_port_rbc.pdf", width = 24, height = 21.6, onefile = T)
grid.arrange(arrangeGrob(plot_townsend_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s70_townsend_port_other.pdf", width = 24, height = 24, onefile = T)
grid.arrange(arrangeGrob(plot_townsend_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()


#================================= GWAS w/ 300K WB =================================
# Read in group level and individual PGS files for covariates
group_pgs_df_300K <- read_tsv("data/pgs_pred/group_non_pgs_df_300K.tsv") %>%
  filter(threshold == 1)
ind_pgs_df_300K <- read_tsv("data/pgs_pred/ind_pgs_df_300K.tsv") %>%
  select(-(16:34))

# Upper limit for the GWAS set
pc_dist_gwas_300K <- read_tsv("data/pca/pc_dist_best_gwas_300K_std.tsv")
upper_300K <- unname(quantile(pc_dist_gwas_300K$pc_dist, 0.975))

# Get positions of the knots by density
temp_300K <- ind_pgs_df_300K %>% filter(phenotype == "Height" & pc_dist > upper_300K)
temp_300K <- temp_300K %>% arrange(pc_dist)
knots_300K <- c(upper_300K,
              temp_300K$pc_dist[round(nrow(temp_300K) / 9)], 
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 2)], 
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 3)], 
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 4)],
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 5)], 
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 6)], 
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 7)], 
              temp_300K$pc_dist[round(nrow(temp_300K) / 9 * 8)])

# Group level plot
plot_group_level_300K <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 38) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots_300K), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(upper, stop, by = step))
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
  }  else if(trait %in% c("Body_Fat_Perc", "Cystatin_C")){
    ymax <- 3.2
  } else if(trait %in% c("BMI", "Eosinophil")){
    ymax <- 3.4
  } else if(trait %in% c("RBC")){
    ymax <- 2.8
  } else if(trait %in% c("Lymphocyte")){
    ymax <- 5.4
  } else if(trait %in% c("WBC", "Weight")){
    ymax <- 3.8
  } else if(trait == "Platelet", "Monocyte", "MCH", "Triglycerides", "LDL"){
    ymax <- 2.6
  } else if(trait %in% c("MCV")){
    ymax <- 2.4
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
      annotate("text", label = "↑", x = 40, y = 1.694, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 82, y = 1.694, size = 8, hjust = 1, 
               vjust = 0.8, family = "Helvetica") +
      annotate(geom = "segment", x = 100, y = 1.05, xend = 107, yend = 1.05, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 279\nindividuals", x = 127, y = 1.15, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Individual level plot
plot_ind_level_300K <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait) %>%
    drop_na(pc_dist, pred_error)
  
  # Plot prediction error relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 38) %>%
    dplyr::summarise(ref_pred_error = mean(pred_error))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = pred_error / ref_pred_error) %>%
    select(-ref_pred_error) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(pc_dist > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(pc_dist, knots = knots_300K), data = plot_df)

  start <- upper
  stop <- max(plot_df$pc_dist)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(pc_dist = seq(start, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = pc_dist, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(size=2.5, alpha=0.1, shape = 16, color = "#beb0d7")+
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
  
  if(trait %in% c("MCV", "MCH", "Eosinophil", "Cystatin_C")){
    ymin <- 0.4
    ymax <- 4.2
  } else if(trait %in% c("Lymphocyte", "WBC")){
    ymin <- 0.18
    ymax <- 2.2
  } else if(trait == "Monocyte"){
    ymin <- 0.35
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
      annotate("text", label = "↑", x = 40, y = 0.45, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 82, y = 0.45, size = 8, hjust = 1, vjust = 0.8,
               family = "Helvetica")
  }
  return(plot)
}

# Make the plots
plot_group_height_300K <- plot_group_level_300K(group_pgs_df_300K, "Height", upper_300K)
plot_ind_height_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Height", upper_300K)

plot_group_cystatin_c_300K <- plot_group_level_300K(group_pgs_df_300K, "Cystatin_C", upper_300K)
plot_ind_cystatin_c_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Cystatin_C", upper_300K)

plot_group_platelet_300K <- plot_group_level_300K(group_pgs_df_300K, "Platelet", upper_300K)
plot_ind_platelet_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Platelet", upper_300K)

plot_group_mcv_300K <- plot_group_level_300K(group_pgs_df_300K, "MCV", upper_300K)
plot_ind_mcv_300K <- plot_ind_level_300K(ind_pgs_df_300K, "MCV", upper_300K)

plot_group_weight_300K <- plot_group_level_300K(group_pgs_df_300K, "Weight", upper_300K)
plot_ind_weight_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Weight", upper_300K)

plot_group_mch_300K <- plot_group_level_300K(group_pgs_df_300K, "MCH", upper_300K)
plot_ind_mch_300K <- plot_ind_level_300K(ind_pgs_df_300K, "MCH", upper_300K)

plot_group_bmi_300K <- plot_group_level_300K(group_pgs_df_300K, "BMI", upper_300K)
plot_ind_bmi_300K <- plot_ind_level_300K(ind_pgs_df_300K, "BMI", upper_300K)

plot_group_rbc_300K <- plot_group_level_300K(group_pgs_df_300K, "RBC", upper_300K)
plot_ind_rbc_300K <- plot_ind_level_300K(ind_pgs_df_300K, "RBC", upper_300K)

plot_group_body_fat_perc_300K <- plot_group_level_300K(group_pgs_df_300K, "Body_Fat_Perc", upper_300K)
plot_ind_body_fat_perc_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Body_Fat_Perc", upper_300K)

plot_group_monocyte_300K <- plot_group_level_300K(group_pgs_df_300K, "Monocyte", upper_300K)
plot_ind_monocyte_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Monocyte", upper_300K)

plot_group_triglycerides_300K <- plot_group_level_300K(group_pgs_df_300K, "Triglycerides", upper_300K)
plot_ind_triglycerides_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Triglycerides", upper_300K)

plot_group_lymphocyte_300K <- plot_group_level_300K(group_pgs_df_300K, "Lymphocyte", upper_300K)
plot_ind_lymphocyte_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Lymphocyte", upper_300K)

plot_group_wbc_300K <- plot_group_level_300K(group_pgs_df_300K, "WBC", upper_300K)
plot_ind_wbc_300K <- plot_ind_level_300K(ind_pgs_df_300K, "WBC", upper_300K)

plot_group_eosinophil_300K <- plot_group_level_300K(group_pgs_df_300K, "Eosinophil", upper_300K)
plot_ind_eosinophil_300K <- plot_ind_level_300K(ind_pgs_df_300K, "Eosinophil", upper_300K)

plot_group_ldl_300K <- plot_group_level_300K(group_pgs_df_300K, "LDL", upper_300K)
plot_ind_ldl_300K <- plot_ind_level_300K(ind_pgs_df_300K, "LDL", upper_300K)

plot_300K_1 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_height_300K, plot_group_height, plot_ind_height_300K, plot_ind_height, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_weight_300K, plot_group_weight, plot_ind_weight_300K, plot_ind_weight,
                          NULL, NULL, NULL, NULL, 
                          plot_group_bmi_300K, plot_group_bmi, plot_ind_bmi_300K, plot_ind_bmi, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_body_fat_perc_300K, plot_group_body_fat_perc, plot_ind_body_fat_perc_300K, plot_ind_body_fat_perc, 
                   labels = c('A. Height (group level, both WB and NWB in the prediction set)\n',
                              'B. Height (group level, only NWB in the prediction set, main\ntext ver.)',
                              'C. Height (individual level, both WB and NWB in the prediction\nset)',
                              'D. Height (individual level, only NWB in the prediction set,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Weight (group level, both WB and NWB in the prediction set)\n',
                              'F. Weight (group level, only NWB in the prediction set, main\ntext ver.)',
                              'G. Weight (individual level, both WB and NWB in the prediction\nset)',
                              'H. Weight (individual level, only NWB in the prediction set,\nmain text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. BMI (group level, both WB and NWB in the prediction set)\n', 
                              'J. BMI (group level, only NWB in the prediction set, main text\nver.)',
                              'K. BMI (individual level, both WB and NWB in the prediction\nset)',
                              'L. BMI (individual level, only NWB in the prediction set, main\ntext ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Body fat percentage (group level, both WB and NWB in the\nprediction set)', 
                              'N. Body fat percentage (group level, only NWB in the\nprediction set, main text ver.)', 
                              'O. Body fat percentage (individual level, both WB and NWB in\nthe prediction set)',
                              'P. Body fat percentage (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_300K_2 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_monocyte_300K, plot_group_monocyte, plot_ind_monocyte_300K, plot_ind_monocyte, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_lymphocyte_300K, plot_group_lymphocyte, plot_ind_lymphocyte_300K, plot_ind_lymphocyte, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_wbc_300K, plot_group_wbc, plot_ind_wbc_300K, plot_ind_wbc, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_eosinophil_300K, plot_group_eosinophil, plot_ind_eosinophil_300K, plot_ind_eosinophil, 
                   labels = c('A. Monocyte count (group level, both WB and NWB in the\nprediction set)', 
                              'B. Monocyte count (group level, only NWB in the prediction\nset, main text ver.)',
                              'C. Monocyte count (individual level, both WB and NWB in the\nprediction set)',
                              'D. Monocyte count (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Lymphocyte count (group level, both WB and NWB in the\nprediction set)', 
                              'F. Lymphocyte count (group level, only NWB in the prediction\nset, main text ver.)',
                              'G. Lymphocyte count (individual level, both WB and NWB in\nthe prediction set)', 
                              'H. Lymphocyte count (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. White blood cell count (group level, both WB and NWB in\nthe prediction set)', 
                              'J. White blood cell count (group level, only NWB in the\nprediction set, main text ver.)',
                              'K. White blood cell count (individual level, both WB and NWB\nin the prediction set)', 
                              'L. White blood cell count (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. Eosinophil count (group level, both WB and NWB in the\nprediction set)', 
                              'N. Eosinophil count (group level, only NWB in the prediction\nset, main text ver.)',
                              'O. Eosinophil count (individual level, both WB and NWB in the\nprediction set)', 
                              'P. Eosinophil count (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_300K_3 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_mcv_300K, plot_group_mcv, plot_ind_mcv_300K, plot_ind_mcv, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_mch_300K, plot_group_mch, plot_ind_mch_300K, plot_ind_mch,
                          NULL, NULL, NULL, NULL, 
                          plot_group_rbc_300K, plot_group_rbc, plot_ind_rbc_300K, plot_ind_rbc, 
                   labels = c('A. Mean corpuscular volume (group level, both WB and NWB in\nthe prediction set)', 
                              'B. Mean corpuscular volume (group level, only NWB in the\nprediction set, main text ver.)',
                              'C. Mean corpuscular volume (individual level, both WB and\nNWB in the prediction set)',
                              'D. Mean corpuscular volume (individual level, only NWB in\nthe prediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Mean corpuscular hemoglobin (group level, both WB and\nNWB in the prediction set)',
                              'F. Mean corpuscular hemoglobin (group level, only NWB in\nthe prediction set, main text ver.)',
                              'G. Mean corpuscular hemoglobin (individual level, both WB\nand NWB in the prediction set)',
                              'H. Mean corpuscular hemoglobin (individual level, only NWB\nin the prediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Red blood cell count (group level, both WB and NWB in the\nprediction set)', 
                              'J. Red blood cell count (group level, only NWB in the\nprediction set, main text ver.)',
                              'K. Red blood cell count (individual level, both WB and NWB in\nthe prediction set)', 
                              'L. Red blood cell count (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 6,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

plot_300K_4 <- plot_grid(NULL, NULL, NULL, NULL, 
                          plot_group_cystatin_c_300K, plot_group_cystatin_c, plot_ind_cystatin_c_300K, plot_ind_cystatin_c, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_platelet_300K, plot_group_platelet, plot_ind_platelet_300K, plot_ind_platelet, 
                          NULL, NULL, NULL, NULL, 
                          plot_group_triglycerides_300K, plot_group_triglycerides, plot_ind_triglycerides_300K, plot_ind_triglycerides,
                          NULL, NULL, NULL, NULL, 
                          plot_group_ldl_300K, plot_group_ldl, plot_ind_ldl_300K, plot_ind_ldl, 
                   labels = c('A. Cystatin C level (group level, both WB and NWB in the\nprediction set)', 
                              'B. Cystatin C level (group level, only NWB in the prediction\nset, main text ver.)',
                              'C. Cystatin C level (individual level, both WB and NWB in the\nprediction set)', 
                              'D. Cystatin C level (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'E. Platelet count (group level, both WB and NWB in the\nprediction set)', 
                              'F. Platelet count (group level, only NWB in the prediction set,\nmain text ver.)', 
                              'G. Platelet count (individual level, both WB and NWB in the\nprediction set)',
                              'H. Platelet count (individual level, only NWB in the prediction\nset, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'I. Triglyceride level (group level, both WB and NWB in the\nprediction set)', 
                              'J. Triglyceride level (group level, only NWB in the prediction\nset, main text ver.)', 
                              'K. Triglyceride level (individual level, both WB and NWB in the\nprediction set)', 
                              'L. Triglyceride level (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              '',
                              'M. LDL cholesterol level (group level, both WB and NWB in the\nprediction set)', 
                              'N. LDL cholesterol level (group level, only NWB in the\nprediction set, main text ver.)',
                              'O. LDL cholesterol level (individual level, both WB and NWB in\nthe prediction set)', 
                              'P. LDL cholesterol level (individual level, only NWB in the\nprediction set, main text ver.)',
                              '',
                              '',
                              '',
                              ''), 
                   ncol = 4, nrow = 8,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.25, 1, 0.25, 1, 0.25, 1, 0.25, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s72_300K_gwas_physical.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_300K_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s73_300K_gwas_wbc.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_300K_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s74_300K_gwas_rbc.pdf", width = 48, height = 20.45, onefile = T)
grid.arrange(arrangeGrob(plot_300K_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s75_300K_gwas_other.pdf", width = 48, height = 27.27, onefile = T)
grid.arrange(arrangeGrob(plot_300K_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

           
plot_group_level_300K_zoomed <- function(pgs_df, trait, upper){
  
  plot_df <- pgs_df %>% filter(phenotype == trait)
  
  # Plot prediction accuracy relative to the 25 bins with a genetic distance most similar
  # to the GWAS set
  denominator <- plot_df %>%
    group_by(phenotype) %>%
    filter(group_close_to_gwas <= 38) %>%
    dplyr::summarise(ref_partial_r2 = mean(partial))
  
  plot_df <- plot_df %>%
    left_join(denominator, by = c("phenotype")) %>%
    mutate(relative_performance = partial / ref_partial_r2) %>%
    select(-ref_partial_r2) %>%
    ungroup()
  
  # Only plot the data points > upper
  plot_df <- plot_df %>% filter(median_pc > upper)
  
  # Fit a spline
  lm <- lm(relative_performance ~ 
             splines::bs(median_pc, knots = knots_300K), data = plot_df)

  start <- upper
  stop <- max(plot_df$median_pc)
  step <- (stop - start) / 400
           
  plot_df_2 <- cbind.data.frame(median_pc = seq(upper, stop, by = step))
  plot_df_2$fitted_val = unname(predict(lm, newdata = plot_df_2))
  
  plot <- ggplot(plot_df, aes(x = median_pc, y = relative_performance)) +
    geom_hline(yintercept = 1, size = 0.5) +
    geom_point(data = subset(plot_df, median_values <= 5.15),
               size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
    scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
    geom_line(data = plot_df_2,
              aes(x = median_pc, y = fitted_val), size = 2.5, color = "#e66100")
  
  # Set the upper limit of the y-axis
  if(trait == "Height"){
    ymax <- 1.8
  }  else if(trait %in% c("Body_Fat_Perc")){
    ymax <- 3.2
  } else if(trait %in% c("BMI", "Eosinophil")){
    ymax <- 3.4
  } else if(trait %in% c("RBC", "Cystatin_C")){
    ymax <- 2.8
  } else if(trait %in% c("Lymphocyte")){
    ymax <- 5.4
  } else if(trait %in% c("WBC", "Weight")){
    ymax <- 3.8
  } else if(trait == "Platelet", "Monocyte", "MCH", "Triglycerides", "LDL"){
    ymax <- 2.6
  } else if(trait %in% c("MCV")){
    ymax <- 2.4
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
    coord_cartesian(xlim = c(upper, 5.2), ylim = c(0.0, ymax))
  
  # Add some labels for height
  if(trait == "Height"){
    plot <- plot +
      annotate("text", label = "↑", x = 0.6, y = 1.694, size = 15, hjust = 1,
               family = "Helvetica") +
      annotate("text", label = "Better prediction", x = 1.6, y = 1.694, size = 8, hjust = 1, 
               vjust = 0.8, family = "Helvetica") +
      annotate(geom = "segment", x = 2.7, y = 1.6, xend = 3, yend = 1.6, 
               color = "black", size = 1) +
      annotate("text", label = "A bin of 279\nindividuals", x = 3.3, y = 1.6, size = 8, 
               family = "Helvetica")
  }
  return(plot)
}

# Make the plots
plot_group_height_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Height", upper_300K)

plot_group_cystatin_c_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Cystatin_C", upper_300K)

plot_group_platelet_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Platelet", upper_300K)

plot_group_mcv_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "MCV", upper_300K)

plot_group_weight_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Weight", upper_300K)

plot_group_mch_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "MCH", upper_300K)

plot_group_bmi_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "BMI", upper_300K)

plot_group_rbc_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "RBC", upper_300K)

plot_group_body_fat_perc_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Body_Fat_Perc", upper_300K)

plot_group_monocyte_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Monocyte", upper_300K)

plot_group_triglycerides_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Triglycerides", upper_300K)

plot_group_lymphocyte_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Lymphocyte", upper_300K)

plot_group_wbc_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "WBC", upper_300K)

plot_group_eosinophil_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "Eosinophil", upper_300K)

plot_group_ldl_300K_zoomed <- plot_group_level_300K_zoomed(group_pgs_df_300K, "LDL", upper_300K)

plot_300K_zoomed_1 <- plot_grid(NULL, NULL, 
                          plot_group_height_300K_zoomed, plot_group_weight_300K_zoomed, 
                          NULL, NULL, 
                          plot_group_bmi_300K_zoomed, plot_group_body_fat_perc_300K_zoomed, 
                   labels = c('A. Height',
                              'B. Weight',
                              '',
                              '',
                              'C. BMI',
                              'D. Body fat percentage',
                              '',
                              ''), 
                   ncol = 2, nrow = 4,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

plot_300K_zoomed_2 <- plot_grid(NULL, NULL, 
                          plot_group_monocyte_300K_zoomed, plot_group_lymphocyte_300K_zoomed, 
                          NULL, NULL, 
                          plot_group_wbc_300K_zoomed, plot_group_eosinophil_300K_zoomed, 
                   labels = c('A. Monocyte count',
                              'B. Lymphocyte count',
                              '',
                              '',
                              'C. White blood cell count',
                              'D. Eosinophil count',
                              '',
                              ''), 
                   ncol = 2, nrow = 4,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

plot_300K_zoomed_3 <- plot_grid(NULL, NULL, 
                          plot_group_mcv_300K_zoomed, plot_group_mch_300K_zoomed, 
                          NULL, NULL, 
                          plot_group_rbc_300K_zoomed, NULL, 
                   labels = c('A. Mean corpuscular volume',
                              'B. Mean corpuscular hemoglobin',
                              '',
                              '',
                              'C. Red blood cell count',
                              '',
                              '',
                              ''), 
                   ncol = 2, nrow = 4,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

plot_300K_zoomed_4 <- plot_grid(NULL, NULL, 
                          plot_group_cystatin_c_300K_zoomed, plot_group_platelet_300K_zoomed, 
                          NULL, NULL, 
                          plot_group_triglycerides_300K_zoomed, plot_group_ldl_300K_zoomed, 
                   labels = c('A. Cystatin C level',
                              'B. Platelet count',
                              '',
                              '',
                              'C. Triglyceride level',
                              'D. LDL cholesterol level',
                              '',
                              ''), 
                   ncol = 2, nrow = 4,
                   label_x = 0.01, hjust = 0,
                   label_size = 28, scale = 1,
                   rel_heights = c(0.1, 1, 0.1, 1),
                   label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_s76_300K_gwas_physical_zoomed.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_300K_zoomed_1,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s77_300K_gwas_wbc_zoomed.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_300K_zoomed_2,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s78_300K_gwas_rbc_zoomed.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_300K_zoomed_3,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

grDevices::cairo_pdf("img/fig_s79_300K_gwas_other_zoomed.pdf", width = 24, height = 12, onefile = T)
grid.arrange(arrangeGrob(plot_300K_zoomed_4,
                         bottom = textGrob("Genetic distance from the GWAS sample", 
                                           gp=gpar(fontfamily = "Helvetica", fontsize=24))))
dev.off()

