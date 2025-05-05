# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# import
data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% filter(file != "6meo")
data <- data %>% unite(col = "gnccn", source_gnccn, target_gnccn, sep = "_", remove = FALSE)

# define conserved
temp <- data %>% filter((all_cc_cxc_para_ck > 0.5) & (all_non_ackr_para_ckr > 0.5))  %>%
  select(gnccn) %>% unique()
temp <- temp$gnccn
data <- data %>% mutate(conserved_rin = case_when(
  gnccn %in% temp ~ "yes",
  !(gnccn %in% temp) ~ "no",
))
rm(temp)

# data <- data %>% filter(file == "7f1t") %>% select(source_gnccn, target_gnccn, conserved_rin)

# count pct
data <- data %>% group_by(file) %>% count(conserved_rin) %>% ungroup()
data <- data %>% group_by(file) %>% mutate(total = sum(n)) %>% ungroup()
data <- data %>% mutate(pct = n / total)
data <- data %>% group_by(conserved_rin) %>% mutate(mean = mean(pct)) %>% ungroup() %>%
  select(conserved_rin, mean) %>% unique()

# plot
hsize <- 1
data %>%
  ggplot(aes(x = 1, y = mean, fill = conserved_rin)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()

ggsave(filename = "pct_rins_cons.pdf", 
       plot = last_plot(), path = "F2/output/",
       width = 5,
       height = 3)



