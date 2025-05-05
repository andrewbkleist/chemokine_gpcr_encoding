# Name:     01_get_t2_logo_ck.R
# Updated:  20221109
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

################################################################################

  # load tier2 sequence positions
  ck.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_all.csv") %>% 
    filter(type != "specific") %>%
    select(source_gnccn)
    
  # # top predictive positions
  # ck.pos <- c("class", ck.pos$source_gnccn, 
  #             "B2.1", "B2.4", "b1b2.8", "B1.5", "b1b2.16", "cxb1.5",
  #             "b1b2.4", "NTc.Cm1", "B3.4", "b1b2.9", "b1b2.6", "H.10", "B1.1", "b2b3.4", "NTc.Cm2",
  #             "B1.7")
    
   ck.pos <- c("class", ck.pos$source_gnccn, 
                "b2b3.4")
    
  # load alignment df
  ck <- read_csv("02_ck_seq/data/processed/ALL_para_df.csv")
  ck <- ck[(names(ck) %in% ck.pos)]
  rm(ck.pos)
    
  # make CC/CXC and uniqte columns
  ck.cc <- ck %>% filter(class == "cc")
  ck.cc <- unite(ck.cc[,2:ncol(ck.cc)], col = seq_string,  sep = "")
  ck.cxc <- ck %>% filter(class == "cxc")
  ck.cxc <- unite(ck.cxc[,2:ncol(ck.cxc)], col = seq_string,  sep = "")
    
  # plot
  p.cc <- ggseqlogo(ck.cc,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  p.cxc <- ggseqlogo(ck.cxc,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(p.cc, p.cxc)
  #rm(ck, ck.cc, ck.cxc, p.cc, p.cxc)
  
 