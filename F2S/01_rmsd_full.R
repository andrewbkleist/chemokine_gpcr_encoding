# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

a <- read_csv("10_structure_pairwise/output/pairwise_chemokine_rmsd.csv")
b <- read_csv("10_structure_pairwise/output/pairwise_receptor_rmsd.csv")
ab <- read_csv("10_structure_pairwise/output/pairwise_complex_rmsd.csv")

a <- a %>% group_by(file1, file2) %>% 
  mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  select(-gnccn, -RMSD) %>% unique()
a$type <- c("ck")

b <- b %>% group_by(file1, file2) %>% 
  mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  select(-gnccn, -RMSD) %>% unique()
b$type <- c("ckr")

ab <- ab %>% group_by(file1, file2) %>% 
  mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  select(-gnccn, -RMSD) %>% unique() # %>% 
  # filter(file1 != "7sk3") %>% 
  # filter(file2 != "7sk3")
ab$type <- c("complex")

data <- rbind(a, b, ab)

data %>%
   # filter(file1 != "7sk3") %>% 
   # filter(file2 != "7sk3") %>%
  ggplot(aes(type, mean)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2) +
  theme_minimal()

ggsave(filename = "rmsd_mean_violin.pdf", 
       plot = last_plot(), path = "F2/output/",
       width = 4,
       height = 5)

