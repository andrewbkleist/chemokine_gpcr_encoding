# Name:     03_motif_pct.R
# Updated:  20230525
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### 1: CHEMOKINE #############################################################
  
  # import, reformat
  data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")  
  
  # select only relevant info (motif, protein, class, pct ortho)
  data <- data %>% 
    filter(pct_ortho>= 0.5) %>%
    select(motif, count_super) %>% unique() %>% count(count_super)
  
  
  data$total <- c(974)
  data <- data %>% mutate(pct = n/total)
  
  # graph
  data %>%
    ggplot(aes(x = "", y= pct)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 
  
  ggsave(filename = "pct_motif_uniqueness_ck.pdf", 
         plot = last_plot(), path = "F5S/output/",
         width = 3,
         height = 3)
  
##### 2: RECEPTOR - NTERM #######################################################
  
  # import, reformat
  data <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")  
  
  # select only relevant info (motif, protein, class, pct ortho)
  data <- data %>% 
    filter(pct_ortho>= 0.5) %>%
    select(motif, count_super) %>% unique() %>% count(count_super)
  data$total <- c(934)
  data <- data %>% mutate(pct = n/total)
  
  # graph
  data %>%
    ggplot(aes(x = "", y= pct)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) 
  
  ggsave(filename = "pct_motif_uniqueness_ckr.pdf", 
         plot = last_plot(), path = "F5S/output/",
         width = 3,
         height = 3)
  
