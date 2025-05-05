# Name:     02_get_t2_logo_ckr.R
# Updated:  20221109
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

################################################################################

  # load tier2 sequence positions
  ckr.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_all.csv") %>%
    filter(type != "specific") %>%
    select(target_gnccn)
  ckr.pos$gn <- c("gn")
  ckr.pos <- ckr.pos %>% unite(target_gnccn, c(gn, target_gnccn), sep = "")
  ckr.pos <- c("class", ckr.pos$target_gnccn, "gn45x52", "gn6x65", "gn4x65",
               "gn1x28", "gn1x24", "5x36", "gn2x66", "gn7x27", "gn3x29", "gn7x35")
  
  # load alignment df
  ckr <- read_csv("03_ckr_seq/data/processed/ALL_classa_df.csv")
  ckr <- ckr[(names(ckr) %in% ckr.pos)]
  rm(ckr.pos)
  ckr <- ckr %>% filter(class %in% c("cc", "cxc"))
  
  # make CC/CXC and uniqte columns
  ckr.cc <- ckr %>% filter(class == "cc")
  ckr.cc <- unite(ckr.cc[,2:ncol(ckr.cc)], col = seq_string,  sep = "")
  ckr.cxc <- ckr %>% filter(class == "cxc")
  ckr.cxc <- unite(ckr.cxc[,2:ncol(ckr.cxc)], col = seq_string,  sep = "")
  
  # plot
  p.cc <- ggseqlogo(ckr.cc,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  p.cxc <- ggseqlogo(ckr.cxc,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  gridExtra::grid.arrange(p.cc, p.cxc)
