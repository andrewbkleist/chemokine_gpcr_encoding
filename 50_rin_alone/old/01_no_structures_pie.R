# Name:     01_no_structures_pie.R
# Updated:  20230206
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: PIE CHART UNIQUE STRUCTURAL CONTACTS ##################################

  # import contacts, select interface only, select Xray only
  res <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file != "6meo") %>%
    count(source_gnccn, target_gnccn)
  colnames(res)[3] <- c("no_pdb")
    
  # add descriptive column for number of contacts
  res <- res %>% mutate(no = case_when(
    no_pdb == 1 ~ "one",
    no_pdb == 2 ~ "two",
    no_pdb == 3 ~ "three",
    no_pdb == 4 ~ "four",
    no_pdb == 5 ~ "five",
    no_pdb == 6 ~ "cons",
    no_pdb == 7 ~ "cons",
    no_pdb == 8 ~ "cons",
    no_pdb == 9 ~ "cons",
    no_pdb == 10 ~ "cons",
    no_pdb == 11 ~ "cons"
  ))
  
  # isolate, count number, fractions
  n.unique <- res %>% unique()
  res.cnt <- res %>% unique() %>% count(no)
  res.cnt <- res.cnt %>% mutate(pct = n / nrow(n.unique))
  
  # PLOT HISTO
  res.cnt$no <- factor(res.cnt$no, levels = c("one", "two", "three", "four", "cons"))
  # res.cnt$no <- factor(res.cnt$no, levels = res.cnt$no[order(res.cnt$n)])
  res.cnt %>%
    ggplot(aes(x = "", y= pct, fill = (no))) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 


