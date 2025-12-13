library(tidyr)
library(dplyr)
library(readr)
library(data.table)
library(ggplot2)
library(patchwork)
library(cowplot)

# Read Fst and PC distance files
pc_dist <- read_tsv("data/pca/pc_dist_best_pred_300K_std.tsv")
pc_dist_gwas <- read_tsv("data/pca/pc_dist_best_gwas_300K_std.tsv")

pc_dist_gwas$group <- "GWAS"
pc_dist$group <- "Pred"

best_pc_dist <- rbind.data.frame(pc_dist_gwas, pc_dist)

# Proxy for centroids of the 3 1000 Genomes populations
proxy <- read_tsv("data/kgp_merged/pca_proxy_300K.tsv")

# Upper limit for the GWAS set
upper <- unname(quantile(pc_dist_gwas$pc_dist, 0.975))

# Distribution of genetic distance
plot <- best_pc_dist %>% 
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
  scale_fill_manual(values = c(GWAS = "#E1BE6A", Test = "#40B0A6")) +
  scale_color_manual(values = c(GWAS = "#E1BE6A", Test = "#40B0A6")) +
  coord_trans(xlim = c(0, 200)) +
  scale_x_continuous(breaks=c(0, 20, 60, 100, 140, 180),
                     expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1e+05, 2e+05, 3e+05),
                     labels = c(expression(0~"\u00D7"~10^5), expression(1~"\u00D7"~10^5), 
                                expression(2~"\u00D7"~10^5), expression(3~"\u00D7"~10^5))) +
  annotate(geom = "segment", x = 10, y = 25000, 
           xend = proxy$pc_dist[1], yend = 0, 
           size = 1, arrow = arrow(length = unit(3, "mm"))) +
  annotate("text", label = "CEU", x = 15, y = 35000, size = 8, hjust = 0.5,
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
plot_inset <- best_pc_dist %>% 
  filter(pc_dist <= 3) %>%
  ggplot(aes(x = pc_dist, group = group, fill = group, color = group)) +
  geom_histogram(binwidth = 0.1, alpha = 0.2) +
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
  scale_fill_manual(values = c(GWAS = "#E1BE6A", Test = "#40B0A6")) +
  scale_color_manual(values = c(GWAS = "#E1BE6A", Test = "#40B0A6")) +
  #annotate(geom = "segment", x = 10, y = 10000, xend = 5, yend = 10000, 
  #         color = "#E1BE6A", size = 1) +
  annotate("text", label = "GWAS", x = 2.8, y = 10^7.5, size = 8, hjust = 1,
           family = "Helvetica", color = "#E1BE6A") +
  #annotate(geom = "segment", x = 25, y = 1200, xend = 20, yend = 1200, 
  #         color = "#40B0A6", size = 1) +
  annotate("text", label = "Prediction", x = 2.8, y = 10^6.75, size = 8, hjust = 1,
           family = "Helvetica", color = "#40B0A6") +
  scale_x_continuous(breaks=c(0, 0.5, 1, 1.5, 2, 2.5, 3),
                     expand = c(0, 0)) +
  coord_cartesian(xlim = c(0.25, 3.05), ylim = c(1, 10^8.5), clip = "off") +
  scale_y_continuous(trans='log10', expand = c(0, 0),
                     breaks = c(1, 10^4, 10^8),
                     labels = c(expression(10^0), expression(10^4), expression(10^8))) +
  annotation_logticks(base = 10, sides = "l", outside = T)

# Combine 2 plots for the distribution of genetic distance
plot_dist <- plot + inset_element(plot_inset, 0.05, 0.2, 0.99, 0.99)

grDevices::cairo_pdf("img/fig_s71_gen_dist_300K.pdf", width = 12, height = 6)
print(plot_dist)
dev.off()
