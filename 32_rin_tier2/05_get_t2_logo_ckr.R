# Name:     05_get_t2_logo_ckr.R
# Updated:  20221122
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

################################################################################

  # define tier2 sequence positions
  pos <- c("gn1x24", "gn6x58", "gn6x62")
  
  # load alignment df
  seq <- read_csv("03_ckr_seq/data/processed/ALL_classa_df.csv")
  
  # make CC/CXC and unique columns
  group1 <- seq %>% filter(class == "cc")
  group1 <- group1[(names(group1) %in% pos)]
  group1 <- unite(group1[,1:ncol(group1)], col = seq_string,  sep = "")
  
  group2 <- seq %>% filter(class == "cxc")
  group2 <- group2[(names(group2) %in% pos)]
  group2 <- unite(group2[,1:ncol(group2)], col = seq_string,  sep = "")
  
  
  # plot
  p1 <- ggseqlogo(group1,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p2 <- ggseqlogo(group2,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(p1, p2)
  #rm(seq, cc, cxc, p.cc, p.cxc)
  