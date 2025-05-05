# Name:     03_str_to_unstr_interactions.R
# Updated:  20201201
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### STRUCTURED TO UNSTRUCTURED COMPLEX ANALYSIS ##############################

  # import contacts
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    select(dom1, dom2, file)
  
  # add labels
  rin <- rin %>% mutate(str_ck = case_when(
    dom1 == "NTc" ~ "NTc",
    dom1 != "NTc" ~ "core_c"
  )) %>% mutate(str_ckr = case_when(
    dom2 == "NTr" ~ "NTr",
    dom2 == "ECL2" ~ "ECL2",
    dom2 != "NTr" & dom2 != "ECL2" ~ "core_r"
  ))
  # note that ECL1 and ECL3 are short loops and are wrapped into the core,
  # similar to chemokine loops
  
  # summary stats
  all <- rin %>% count(str_ck, str_ckr, file)
  all <- all %>% group_by(file) %>% mutate(total_per_file = sum(n)) %>% ungroup()
  all <- all %>% mutate(pct = n / total_per_file)
  all <- all %>% select(-file, -n, -total_per_file) %>% 
    group_by(str_ck, str_ckr) %>% summarise(mean = mean(pct)) %>% ungroup()

  # plot
  # plotted manually in illustrator by varying line width according to rounded percentages
        

  