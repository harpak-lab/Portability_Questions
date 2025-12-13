library(tidyverse)

fst_df <- data.frame()

for (file in list.files(path = 'data/fst',
                        pattern = 'fst_[0-9]+.est', full.names = T)) {
  
  fst_df_append <- read_tsv(file, col_names="Fst",col_types="c")
  fst_df_append$Fst <- fst_df_append$Fst %>% parse_number()
  fst_df_append$label <- rep(c("IID", "Mean_Fst", "Weighted_Fst"), nrow(fst_df_append) / 3)
  fst_df <- rbind(fst_df, fst_df_append)
  rm(fst_df_append)
  
}

row_num = c()
for(i in 1:round(nrow(fst_df) / 3)){
  row_num <- c(row_num, rep(i, 3))
}

fst_df <- fst_df %>% mutate(row_num = row_num) %>%
  pivot_wider(names_from = label, values_from = Fst) %>% select(-row_num) %>%
  arrange(IID)

fst_df %>% write_tsv('data/fst/final_fst.tsv')
