library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)

fst_df <- read_tsv('data/fst/final_fst.tsv')
pc <- read_table('data/extracted_data_fields/pc.txt')
cn <- c()
for(i in 1:40){
  cn[i] <- paste0("PC", i)
}

colnames(pc) <- c("IID", cn)

fst_df <- fst_df %>% left_join(pc, by = "IID")

race <- read_table('data/extracted_data_fields/ethnic_background.txt')
colnames(race) <- c("IID", "race_1", "race_2", "race_3")
race <- race %>% filter(IID %in% fst_df$IID)
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

fst_df <- fst_df %>% left_join(race, by = "IID")
fst_df$race <- factor(fst_df$race, 
                      levels = c("White", "Mixed", "Asian", "Black", "Other", "NA"))

corr <- c(cor(fst_df$Weighted_Fst, fst_df$PC1), cor(fst_df$Weighted_Fst, fst_df$PC2),
          cor(fst_df$Weighted_Fst, fst_df$PC3), cor(fst_df$Weighted_Fst, fst_df$PC40))

plot_fst_pc_1 <- fst_df %>%
  ggplot(aes(x = Weighted_Fst, y = PC1, color = race)) + 
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  xlab(expression(F[st]~between~individual~and~GWAS~sample)) +
  ylab("PC 1") +
  ggtitle(paste("R =", round(corr[1], 4))) +
  theme_bw() +
  theme(legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        legend.position=c(0.9, 0.14),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        plot.title=element_text(size=28, family = "Helvetica")) +
  guides(color = guide_legend(title = "Self-reported race", ncol = 2))

plot_fst_pc_2 <- fst_df %>%
  ggplot(aes(x = Weighted_Fst, y = PC2, color = race)) + 
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  xlab(expression(F[st]~between~individual~and~GWAS~sample)) +
  ylab("PC 2") +
  ggtitle(paste("R =", round(corr[2], 4))) +
  theme_bw() +
  theme(legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        legend.position=c(0.9, 0.14),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        plot.title=element_text(size=28, family = "Helvetica")) +
  guides(color = guide_legend(title = "Self-reported race", ncol = 2))

plot_fst_pc_3 <- fst_df %>%
  ggplot(aes(x = Weighted_Fst, y = PC3, color = race)) + 
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  xlab(expression(F[st]~between~individual~and~GWAS~sample)) +
  ylab("PC 3") +
  ggtitle(paste("R =", round(corr[3], 4))) +
  theme_bw() +
  theme(legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        legend.position=c(0.9, 0.14),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        plot.title=element_text(size=28, family = "Helvetica")) +
  guides(color = guide_legend(title = "Self-reported race", ncol = 2))

plot_fst_pc_40 <- fst_df %>%
  ggplot(aes(x = Weighted_Fst, y = PC40, color = race)) + 
  geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "lm", se = F, color = "black") +
  xlab(expression(F[st]~between~individual~and~GWAS~sample)) +
  ylab("PC 40") +
  ggtitle(paste("R =", round(corr[4], 4))) +
  theme_bw() +
  theme(legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        legend.position=c(0.9, 0.14),
        panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        plot.title=element_text(size=28, family = "Helvetica")) +
  guides(color = guide_legend(title = "Self-reported race", ncol = 2))

pdf(file = "img/fig_s30_pc_1_vs_pc_fst.pdf", width = 12, height = 6)
print(plot_fst_pc_1)
dev.off()

pdf(file = "img/fig_s31_pc_2_vs_pc_fst.pdf", width = 12, height = 6)
print(plot_fst_pc_2)
dev.off()

pdf(file = "img/fig_s32_pc_3_vs_pc_fst.pdf", width = 12, height = 6)
print(plot_fst_pc_3)
dev.off()

pdf(file = "img/fig_s33_pc_40_vs_pc_fst.pdf", width = 12, height = 6)
print(plot_fst_pc_40)
dev.off()
