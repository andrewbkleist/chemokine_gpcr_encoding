# Name:     5_functional_data.R
# Updated:  20191030
# User:     Andrew Kleist
# Figure:   Figure 5H,I

##### LOAD PACKAGES, SET WD ####################################################
  
  # packages, working directory
  library(tidyverse)
  require(Biostrings)
  library(UpSetR)
  library(reshape2)
  setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/kleist_2020/5_ckr_motif/")


##### 1: CXCR4 SULFOTYROSINE DATA ##############################################
  
  data <- read_csv("input/ziarek_2013.csv") %>% select(variant, max_rat, ec50_rat)
  
  order <- c("WT", "Y7A", "Y12A", "Y21A")
  data$variant <- factor(data$variant, levels = rev(order))
  
  data %>%
    ggplot(aes(variant, max_rat)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    coord_flip()
  
  data %>%
    ggplot(aes(variant, ec50_rat)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    coord_flip()
  
##### 2: CCR8 ECL2 MUTATION DATA ###############################################
  
  data <- read_csv("input/barington_2016.csv")
  data %>%
    ggplot(aes(mut, fold_dec_ec50_wt)) +
    geom_bar(stat="identity") +
    # ylim(0,15) +
    theme_minimal()
  
  
    
    