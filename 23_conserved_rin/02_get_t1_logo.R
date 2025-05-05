# Name:     07_get_t1_logo.R
# Updated:  20210330
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

##### 1: TIER 1 SEQUENCE LOGOS #################################################

  # (1.1) CHEMOKINE ------------------------------------------------------------

  # manually define positions for seq logo
  pos <- c("B3.3")
  
  # load alignment df
  ck <- read_csv("02_ck_seq/data/processed/ALL_para_df.csv")
  ck <- ck[(names(ck) %in% pos)]
  ck <- unite(ck[,1:ncol(ck)], col = seq_string,  sep = "")

  p <- ggseqlogo(ck,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(p)

  
  # (1.2) RECEPTOR -------------------------------------------------------------
  
  # manually define positions for seq logo
  pos <- c("gnNTr.Cm1")
  
  # load alignment df
  ck <- read_csv("03_ckr_seq/data/processed/ALL_non_ackr_para_df.csv")
  ck <- ck[(names(ck) %in% pos)]
  ck <- unite(ck[,1:ncol(ck)], col = seq_string,  sep = "")
  
  p <- ggseqlogo(ck,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(p)
  
  