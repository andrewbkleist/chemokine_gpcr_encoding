# Name:     02_custon_motif_heatmaps.R
# Updated:  20210519
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### 1: CUSTOM HEATMAPS #######################################################
  
  # (1) SULFOTYROSINE MOTIFS
  # import, reformat
  data <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")
  
  # select conserved
  # data <- data %>% filter(pct_ortho >= 0.5)
  
  # select only relevant info (motif, protein, class, pct ortho)
  data <- data %>% select(motif, protein, class, pct_ortho)
  
  # order receptors
  order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                   "ccr5","ccr6","ccr7","ccr8","ccr9",
                                   "ccr10","cxcr1","cxcr2","cxcr3",
                                   "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                   "ackr1","ackr2","ackr3","ackr4",
                                   "ccrl2")))
  levels(data$protein)
  data$protein <- factor(data$protein, levels = order.ckr)
  
  # select MOTIF "SLICE"
  data <- data %>% filter(motif == "DY" | motif == "YD" | 
                            motif == "EY" | motif == "YE" | 
                            motif == "DxY" | motif == "DxxY" | 
                            motif == "YxxD" | motif == "YxD" |
                            motif == "ExY" | motif == "ExxY" | 
                            motif == "YxxE" | motif == "YxE")
  
  order.mot <- as.factor(unique(c("DY","YD", "EY", "YE",
                                  "DxY","YxD", "ExY", "YxE",
                                  "DxxY", "YxxD", "ExxY", "YxxE")))
  data$motif <- factor(data$motif, levels = rev(order.mot))
  
  data %>%
    ggplot() + 
    geom_tile(aes(protein, motif), fill = "mediumorchid4") +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
