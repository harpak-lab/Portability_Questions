library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)

# Read Fst and PC distance files
fst <- read_tsv('data/fst/final_fst.tsv')
pc_dist <- read_tsv("data/pca/pc_dist_best_pred_std.tsv")
pc_dist_further_rand_id <- read_tsv("data/ukb_populations/pc_futher_fst_id.txt",
                                    col_names = "IID")
pc_dist_further_rand_id <- pc_dist_further_rand_id$IID
pc_dist_further_rand <- pc_dist %>% filter(IID %in% pc_dist_further_rand_id)
pc_dist_further_rand <- pc_dist_further_rand %>% left_join(fst, by = "IID")

# Convert race codings to labels
race <- read_table("data/extracted_data_fields/ethnic_background.txt")
colnames(race) <- c("IID", "race_1", "race_2", "race_3")
race = race %>% filter(IID %in% pc_dist_further_rand_id)
race$race_1 <- ifelse(race$race_1 < 0, NA, race$race_1)
race$race_1 <- ifelse(is.na(race$race_1), race$race_2, race$race_1)
race$race_1 <- ifelse(race$race_1 < 0, NA, race$race_1)
race$race_1 <- ifelse(is.na(race$race_1), race$race_3, race$race_1)
race$race_1 <- ifelse(race$race_1 < 0, NA, race$race_1)
race <- race %>% select(IID, race_1)
race$race <- ifelse(race$race_1 %in% c(1, 1001, 1002, 1003), "White",
                   ifelse(race$race_1 %in% c(2, 2001, 2002, 2003, 2004), "Mixed",
                          ifelse(race$race_1 %in% c(3, 3001, 3002, 3003, 3004, 5), "Asian",
                                 ifelse(race$race_1 %in% c(4, 4001, 4002, 4003), "Black",
                                        ifelse(race$race_1 == 6, "Other", "NA")))))

race$race <- ifelse(is.na(race$race_1), "NA", race$race)
race <- race %>% select(-race_1)

# Join PC distance and Fst
pc_dist_further_rand <- pc_dist_further_rand %>% left_join(race, by = "IID")
pc_dist_further_rand$race <- factor(pc_dist_further_rand$race, 
                                   levels = c("White", "Mixed", "Asian", "Black", "Other", "NA"))

# Correlate PC distance and Fst
plot_pc_fst <- pc_dist_further_rand %>% 
  ggplot(aes(x = pc_dist, y = Weighted_Fst)) + 
  geom_point(alpha = 0.3, size = 3, aes(color = race)) +
  xlab("Genetic distance from the GWAS sample") +
  ylab(expression(atop(F[st]~between~individual, and~GWAS~sample))) +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        legend.position=c(0.92, 0.34),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_smooth(se = F, method = "lm", color = "black") +
  scale_shape_manual(values=seq(0,19)) +
  coord_trans(xlim = c(0, 200)) +
  scale_x_continuous(breaks=c(0, 20, 60, 100, 140, 180),
                     expand = c(0, 0)) +
  guides(color = guide_legend(title = "Self-reported\nrace"))

# Combined GWAS and prediction sets
pc_dist_gwas <- read_tsv("data/pca/pc_dist_best_gwas_std.tsv")

pc_dist_gwas$group <- "GWAS"
pc_dist$group <- "Pred"

best_pc_dist <- rbind.data.frame(pc_dist_gwas, pc_dist)

# Proxy for centroids of the 3 1000 Genomes populations
proxy <- read_tsv("data/kgp_merged/pca_proxy.tsv")

# Upper limit for the GWAS set
upper <- unname(quantile(pc_dist_gwas$pc_dist, 0.975))

# Distribution of genetic distance
plot_dist <- best_pc_dist %>% 
  ggplot(aes(x = pc_dist, group = group, fill = group, color = group)) +
  geom_histogram(binwidth = 1, alpha = 0.2) +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = upper, size = 1, linetype = "dashed") +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Count") +
  guides(fill = "none",
         color = "none") +
  scale_fill_manual(values = c(GWAS = "#E1BE6A", Pred = "#40B0A6")) +
  scale_color_manual(values = c(GWAS = "#E1BE6A", Pred = "#40B0A6")) +
  coord_trans(xlim = c(0, 200)) +
  scale_x_continuous(breaks=c(0, 20, 60, 100, 140, 180),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1e+05, 2e+05, 3e+05),
                     labels = c(expression(0~"\u00D7"~10^5), expression(1~"\u00D7"~10^5), 
                                expression(2~"\u00D7"~10^5), expression(3~"\u00D7"~10^5))) +
  annotate(geom = "segment", x = 15, y = 25000, 
           xend = proxy$pc_dist[1], yend = 0, 
           size = 1, arrow = arrow(length = unit(3, "mm"))) +
  annotate("text", label = "CEU", x = 22, y = 35000, size = 8, hjust = 0.5,
           family = "Helvetica") +
  annotate(geom = "segment", x = proxy$pc_dist[2], y = 20000, 
           xend = proxy$pc_dist[2], yend = 0, 
           size = 1, arrow = arrow(length = unit(3, "mm"))) +
  annotate("text", label = "CHB", x = proxy$pc_dist[2], y = 35000, size = 8, hjust = 0.5,
           family = "Helvetica") +
  annotate(geom = "segment", x = proxy$pc_dist[3], y = 20000, 
           xend = proxy$pc_dist[3], yend = 0, 
           size = 1, arrow = arrow(length = unit(3, "mm"))) +
  annotate("text", label = "YRI", x = proxy$pc_dist[3], y = 35000, size = 8, hjust = 0.5,
           family = "Helvetica")

# Inset for distribution of genetic distance
plot_dist_inset <- best_pc_dist %>% 
  ggplot(aes(x = pc_dist, group = group, fill = group, color = group)) +
  geom_histogram(binwidth = 1, alpha = 0.2) +
  theme_light() + 
  theme(axis.title=element_blank(),
        axis.text=element_text(size=20, family = "Helvetica"),
        axis.text.y = element_text(vjust = 0.25, hjust = 0.5),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_vline(xintercept = upper, size = 1, linetype = "dashed") +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Count\n") +
  guides(fill = "none",
         color = "none") +
  scale_fill_manual(values = c(GWAS = "#E1BE6A", Pred = "#40B0A6")) +
  scale_color_manual(values = c(GWAS = "#E1BE6A", Pred = "#40B0A6")) +
  annotate(geom = "segment", x = 10, y = 10000, xend = 5, yend = 10000, 
           color = "#E1BE6A", size = 1) +
  annotate("text", label = "GWAS", x = 20, y = 10000, size = 8, hjust = 1,
           family = "Helvetica", color = "#E1BE6A") +
  annotate(geom = "segment", x = 25, y = 1200, xend = 20, yend = 1200, 
           color = "#40B0A6", size = 1) +
  annotate("text", label = "Prediction", x = 40, y = 1200, size = 8, hjust = 1,
           family = "Helvetica", color = "#40B0A6") +
  scale_x_continuous(breaks=c(0, 20, 60, 100, 140, 180),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(0, 200), ylim = c(1, 22000)) +
  scale_y_continuous(trans='log10', expand = c(0, 0),
                     breaks = c(1, 100, 10000),
                     labels = c(expression(10^0), expression(10^2), expression(10^4))) +
  annotation_logticks(base = 10, sides = "l")

# Combine 2 plots for the distribution of genetic distance
plot_dist <- plot_dist + inset_element(plot_dist_inset, 0.05, 0.2, 0.99, 0.99)

# Combine all the panels
# Leave panel A blank and later edit it in Adobe Illustrator
fig_1_a_b <- plot_grid(NULL, NULL, NULL, plot_pc_fst,
                      labels = c('A.', 
                                 'B.',
                                 '',
                                 ''), ncol = 2, nrow = 2,
                      label_x = 0.01, hjust = 0,
                      label_size = 28, scale = 1,
                      rel_heights = c(0.1, 1),
                      label_fontfamily = "Helvetica")

fig_1 <- plot_grid(fig_1_a_b, NULL, plot_dist,
                  labels = c('', 
                             'C.',
                             ''), ncol = 1, nrow = 3,
                  label_x = 0.01, hjust = 0,
                  label_size = 28, scale = 1,
                  rel_heights = c(1, 0.1, 1.5),
                  label_fontfamily = "Helvetica")

grDevices::cairo_pdf("img/fig_1_method.pdf", width = 24, height = 15)
print(fig_1)
dev.off()
