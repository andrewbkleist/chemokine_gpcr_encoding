# Name:     01_tier1_scores.R
# Updated:  20201108
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

##### 1: FETCH T1 SCORES CK AND CKR ############################################

  # CHEMOKINE
  ck <- read_csv("02_ck_seq/output/CK_CONSERVATION.csv") %>% select(ccn, all_ortho) %>% unique()
  ck$ccn <- factor(ck$ccn, levels = ck$ccn[order(ck$all_ortho)])
  
  # plot top
  ck %>%
    top_n(10, all_ortho) %>%
    #filter(all_ortho >= .8) %>%
    ggplot(aes(ccn, all_ortho)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  
  # RECEPTOR
  ckr <- read_csv("03_ckr_seq/output/CKR_CONSERVATION.csv") %>% 
    separate(gn, into = c("temp", "gn"), sep = "gn") %>% select(-temp) %>% 
    select(gn, all_ortho) %>% unique()
  ckr$gn <- factor(ckr$gn, levels = ckr$gn[order(ckr$all_ortho)])
  
  # subset by interface in chemokine-GPCR complexes
  interface <- read_csv("01_structure_contacts/output/RIN_residue.csv") %>% filter(Chain1 != Chain2) %>%
    select(target_gnccn) %>% unique()
  
  interface <- unique(interface$target_gnccn)
  
  # plot top
  ckr %>%
    filter(gn %in% interface) %>%
    top_n(10, all_ortho) %>%
    #filter(mean >= .75) %>%
    ggplot(aes(gn, all_ortho)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  
  