# Name:     01_tier2_scores_rank_boxplot.R
# Updated:  20221109
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: FETCH T2 SCORES CK AND CKR ############################################

  # CHEMOKINE
  ck <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
  ck$motif <- factor(ck$motif, levels = ck$motif[order(ck$mean)])
  
  # plot top
  ck %>%
    #top_n(20, mean) %>%
    filter(mean >= .8) %>%
    ggplot(aes(motif, mean)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  
  # RECEPTOR
  ckr <- read_csv("03_ckr_seq/output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>% 
    separate(motif, into = c("temp", "motif"), sep = "gn") %>% select(-temp)
  ckr$motif <- factor(ckr$motif, levels = ckr$motif[order(ckr$mean)])
  
  # subset by interface in chemokine-GPCR complexes
  interface <- read_csv("01_structure_contacts/output/RIN_residue.csv") %>% filter(Chain1 != Chain2) %>%
    select(target_gnccn) %>% unique()
  
  interface <- unique(interface$target_gnccn)
  
  # plot top
  ckr %>%
    filter(motif %in% interface) %>%
    #top_n(20, mean) %>%
    filter(mean >= .75) %>%
    ggplot(aes(motif, mean)) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(0.05)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
    coord_flip() +
    ylim(.49,1) +
    theme_minimal()
  