# Updated:  20230123
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
#   filter(file != "6meo") %>%
#   select(source_gnccn, target_gnccn, dom1, dom2) %>%
#   unique()

data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(no_pdb == 1) %>%
  filter(all_para_ck < 0.5) %>%
  filter(all_para_ckr < 0.5) %>%
  filter(file != "6meo") %>%
  filter(file != "ngo") %>%
  
  select(source_gnccn, target_gnccn, dom1, dom2) %>%
  unique() # 176 contacts (of 478 unique contacts = 37%)

data <- data %>%
  mutate(str_unstr_ck = case_when(
    dom1 == "NTc" ~ "unstr",
    dom1 != "NTc" ~ "str",
  )) %>% mutate(str_unstr_ckr = case_when(
    dom2 == "NTr" ~ "unstr",
    dom2 != "NTr" ~ "str",
  )) %>% count(str_unstr_ck, str_unstr_ckr) 

data <- data %>%
  unite(col = type, c(str_unstr_ck, str_unstr_ckr), sep = "_")

data <- data %>% mutate(pct = n / sum(n))

data %>%
  ggplot(aes(x = "", y= pct)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

  