library(dplyr)
library(tidyr)
library(readr)
library(ggplot2)
library(data.table)
library(RColorBrewer)

gwas_id <- read_table('data/ukb_populations/wb_gwas_id.txt')
pred_id <- read_table('data/ukb_populations/nwb_all_id.txt')
pred_id <- pred_id %>% arrange(IID)

# Get bin info
non_pgs_df <- read_tsv('data/pgs_pred/group_non_pgs_df.tsv')
non_pgs_df <- non_pgs_df %>%
  select(IID, pc_dist, weighted_pc_groups) %>%
  distinct()

# Age
age <- read_table('data/extracted_data_fields/age.txt')
age_gwas_mean <- mean(age$`21022-0.0`[age$eid %in% gwas_id$IID], na.rm = T)
age_gwas_sd <- sd(age$`21022-0.0`[age$eid %in% gwas_id$IID], na.rm = T)
# 278 is the number of people in each bin
# Use this number so the SE of GWAS set is comparable to other bins
age_gwas_se <- age_gwas_sd / sqrt(278)

age_plot_df <- age %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid, age = `21022-0.0`) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups) %>%
  summarize(median_pc_dist = median(pc_dist, na.rm = T),
            mean_age = mean(age, na.rm = T),
            sd_age = sd(age, na.rm = T),
            n = sum(!is.na(age)),
            se_age = sd_age / sqrt(n))

plot_age <- age_plot_df %>%
  ggplot(aes(x = median_pc_dist, y = mean_age)) +
  geom_vline(xintercept = 1.913934, linetype = "dashed") +
  geom_hline(yintercept = age_gwas_mean) +
  geom_hline(yintercept = age_gwas_mean - age_gwas_se, linetype = "dashed") +
  geom_hline(yintercept = age_gwas_mean + age_gwas_se, linetype = "dashed") +
  geom_point(size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
  geom_errorbar(aes(ymin = mean_age - se_age, ymax = mean_age + se_age), color = "#ff8934") +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Age (mean ± SE)")

# Townsend
townsend <- read_table('data/extracted_data_fields/townsend.txt')
townsend_gwas_mean <- mean(townsend$`189-0.0`[townsend$eid %in% gwas_id$IID], na.rm = T)
townsend_gwas_sd <- sd(townsend$`189-0.0`[townsend$eid %in% gwas_id$IID], na.rm = T)
townsend_gwas_se <- townsend_gwas_sd / sqrt(278)

townsend_plot_df <- townsend %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid, townsend = `189-0.0`) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups) %>%
  summarize(median_pc_dist = median(pc_dist, na.rm = T),
            mean_townsend = mean(townsend, na.rm = T),
            sd_townsend = sd(townsend, na.rm = T),
            n = sum(!is.na(townsend)),
            se_townsend = sd_townsend / sqrt(n))

plot_townsend <- townsend_plot_df %>%
  ggplot(aes(x = median_pc_dist, y = mean_townsend)) +
  geom_vline(xintercept = 1.913934, linetype = "dashed") +
  geom_hline(yintercept = townsend_gwas_mean) +
  geom_hline(yintercept = townsend_gwas_mean - townsend_gwas_se, linetype = "dashed") +
  geom_hline(yintercept = townsend_gwas_mean + townsend_gwas_se, linetype = "dashed") +
  geom_point(size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
  geom_errorbar(aes(ymin = mean_townsend - se_townsend, ymax = mean_townsend + se_townsend), 
                color = "#ff8934") +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Townsend Deprivation Index\n(mean ± SE)")


# Sex
sex <- read_table('data/extracted_data_fields/genetic_sex.txt')
# Count frequency of females
sex_gwas_mean <- length(sex$`22001-0.0`[sex$eid %in% gwas_id$IID & sex$`22001-0.0` == 0]) / 336923

sex_plot_df <- sex %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid, sex = `22001-0.0`) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups) %>%
  summarize(median_pc_dist = median(pc_dist, na.rm = T),
            n = sum(!is.na(sex)),
            mean_sex = sum(sex == 0) / n)

plot_sex <- sex_plot_df %>%
  ggplot(aes(x = median_pc_dist, y = mean_sex)) +
  geom_vline(xintercept = 1.913934, linetype = "dashed") +
  geom_hline(yintercept = sex_gwas_mean) +
  geom_point(size = 5,alpha = 0.4, fill = "#ff8934", color = "#f8766d", shape = 23)+
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Sex (frequency of females)")

# Income
income <- read_table('data/extracted_data_fields/income.txt')
income <- income %>% filter(!is.na(`738-0.0`) & `738-0.0` != -1 & `738-0.0` != -3)
income$income <- ifelse(income$`738-0.0` == 1, "<18K",
                        ifelse(income$`738-0.0` == 2, "18K-31K",
                               ifelse(income$`738-0.0` == 3, "31K-52K",
                                      ifelse(income$`738-0.0` == 4, "52K-100K", ">100K"))))
income$income <- factor(income$income, levels = c("<18K", "18K-31K", "31K-52K", "52K-100K", ">100K"))
income_gwas_mean <- income %>% filter(eid %in% gwas_id$IID) %>%
  group_by(income) %>%
  summarize(mean_income = n() / 301443)

income_plot_df <- income %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid) %>%
  select(IID, income) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups, income) %>%
  summarize(n = sum(!is.na(income))) %>%
  na.omit(income)

income_per_bin <- income %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid) %>%
  select(IID, income) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups) %>%
  summarize(median_pc_dist = median(pc_dist, na.rm = T),
            total_n = sum(!is.na(income)))

income_plot_df <- income_plot_df %>%
  left_join(income_per_bin, by = "weighted_pc_groups")

income_plot_df$mean_income <- income_plot_df$n / income_plot_df$total_n

plot_income <- income_plot_df %>%
  ggplot(aes(x = median_pc_dist, y = mean_income, fill = income)) +
  geom_hline(yintercept = income_gwas_mean$mean_income[income_gwas_mean$income == ">100K"] ,
             color = "#BD0026") +
  geom_hline(yintercept = income_gwas_mean$mean_income[income_gwas_mean$income == ">100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "52K-100K"] ,
             color = "#F03B20") +
  geom_hline(yintercept = income_gwas_mean$mean_income[income_gwas_mean$income == ">100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "52K-100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "31K-52K"] ,
             color = "#FD8D3C") +
  geom_hline(yintercept = income_gwas_mean$mean_income[income_gwas_mean$income == ">100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "52K-100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "31K-52K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "18K-31K"] ,
             color = "#FEB24C") +
  geom_hline(yintercept = income_gwas_mean$mean_income[income_gwas_mean$income == ">100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "52K-100K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "31K-52K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "18K-31K"] +
               income_gwas_mean$mean_income[income_gwas_mean$income == "<18K"] ,
             color = "#FED976") +
  geom_col(position = "stack", width = 0.5) +
  geom_vline(xintercept = 1.913934, linetype = "dashed") +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Income\n(frequency of each category)") +
  labs(fill = "Household\nincome (£)") +
  scale_fill_manual(values = brewer.pal(6, "YlOrRd")[2:6])

# Country
# Use assessment center and convert into country
country <- read_table('data/extracted_data_fields/assessment_center.txt')
country$country <- ifelse(country$`54-0.0` %in% c(11012, 11021, 11011, 11008, 11024, 11020, 11018, 11010, 11016, 11001, 11017, 11009, 11013, 11002, 11007, 11014, 10003, 11006, 11025, 11026, 11027, 11028), "England",
                          ifelse(country$`54-0.0` %in% c(11005, 11004), "Scotland",
                                 ifelse(country$`54-0.0` %in% c(11003, 11022, 11023), "Wales", NA)))

country$country <- factor(country$country, levels = c("England", "Scotland", "Wales"))

country_gwas_mean <- country %>% filter(eid %in% gwas_id$IID) %>%
  group_by(country) %>%
  summarize(mean_country = n() / 336923)

country_plot_df <- country %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid) %>%
  select(IID, country) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups, country) %>%
  summarize(n = sum(!is.na(country)))

country_per_bin <- country %>% 
  filter(eid %in% pred_id$IID) %>%
  rename(IID = eid) %>%
  select(IID, country) %>%
  right_join(non_pgs_df, by = "IID") %>%
  group_by(weighted_pc_groups) %>%
  summarize(median_pc_dist = median(pc_dist, na.rm = T),
            total_n = sum(!is.na(country)))

country_plot_df <- country_plot_df %>%
  left_join(country_per_bin, by = "weighted_pc_groups")

country_plot_df$mean_country <- country_plot_df$n / country_plot_df$total_n

plot_country <- country_plot_df %>%
  ggplot(aes(x = median_pc_dist, y = mean_country, fill = country)) +
  geom_hline(yintercept = country_gwas_mean$mean_country[country_gwas_mean$country == "Wales"],
             color = "#E31A1C") +
  geom_hline(yintercept = country_gwas_mean$mean_country[country_gwas_mean$country == "Wales"] +
               country_gwas_mean$mean_country[country_gwas_mean$country == "Scotland"] ,
             color = "#FD8D3C") +
  geom_hline(yintercept = income_gwas_mean$mean_income[income_gwas_mean$country == "Wales"] +
               country_gwas_mean$mean_country[country_gwas_mean$country == "Scotland"] +
               country_gwas_mean$mean_country[country_gwas_mean$country == "England"] ,
             color = "#FECC5C") +
  geom_col(position = "stack", width = 0.5) +
  geom_vline(xintercept = 1.913934, linetype = "dashed") +
  theme_bw() + 
  theme(axis.title=element_text(size=24, family = "Helvetica"),
        axis.text=element_text(size=20, family = "Helvetica"),
        legend.title=element_text(size=16, family = "Helvetica"),
        legend.text = element_text(size=16, family = "Helvetica"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120, 140, 160, 180), expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  xlab("Genetic distance from the GWAS sample") +
  ylab("Country\n(frequency of each category)") +
  labs(fill = "Country") +
  scale_fill_manual(values = brewer.pal(4, "YlOrRd")[2:4])

pdf(file = "img/fig_s25_townsend.pdf", width = 12, height = 6, onefile = T)
print(plot_townsend)
dev.off()

pdf(file = "img/fig_s26_income.pdf", width = 12, height = 6, onefile = T)
print(plot_income)
dev.off()

pdf(file = "img/fig_s27_sex.pdf", width = 12, height = 6, onefile = T)
print(plot_sex)
dev.off()

pdf(file = "img/fig_s28_age.pdf", width = 12, height = 6, onefile = T)
print(plot_age)
dev.off()

pdf(file = "img/fig_s29_country.pdf", width = 12, height = 6, onefile = T)
print(plot_country)
dev.off()
