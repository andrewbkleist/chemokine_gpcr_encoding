# Name:     05_get_vmip_t2_t2.R
# Updated:  20230109
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: IDENTIFY TIER 2 TO TIER 2 CONTACTS in 4RWS ############################
  
  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file == "4rws")
  
  # plot t2-by-t2
  rin %>%
    ggplot(aes(cc_cxc_lr_ckr, cc_cxc_lr_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
    #scale_size(range = c(2, 8)) +
    xlim(0.5,1) +
    ylim(0.5,1) +
    theme_minimal()  
  
  # plot t2 score-by-t2 score
  rin %>%
    filter(cc_cxc_lr_ckr >= 0.75) %>%
    filter(cc_cxc_lr_ck >= 0.75) %>%
    ggplot(aes(cc_cxc_lr_score_ckr, cc_cxc_lr_score_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
    #scale_size(range = c(2, 8)) +
    xlim(0,1) +
    ylim(0,1) +
    theme_minimal()  
  
  # define t2-by-t2
  test <- rin %>% filter(cc_cxc_lr_ck >= 0.75 & cc_cxc_lr_ckr >= 0.75)
  test <- rin %>% 
    filter(cc_cxc_lr_ck >= 0.75 & cc_cxc_lr_ckr >= 0.75) %>%
    filter(cc_cxc_lr_score_ck >= 0.75 & cc_cxc_lr_score_ckr >= 0.75)
  
  
  # remove distal N-term contacts
  rin.t2 <- rin.t2 %>% filter(source_gnccn != "NTc.Cm10") %>% 
    filter(source_gnccn != "NTc.Cm8") %>% filter(source_gnccn != "NTc.Cm7")
  
  # write CONECT
  WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
                          rin.t2,
                          "4rws",
                          "01_structure_contacts/data/pdbs/4rws_ck_clean.pdb",
                          "35_viral/output/4rws_t2_t2.csv")
  
  
  
##### 2: CC/CXC SEQUENCE LOGOS AT vMIP-II-CXCR4 INTERACTING POSITIONS ##########
  
  # (1.1) CHEMOKINE ------------------------------------------------------------
  
  ck.pos <- c("class", "NTc.Cm1", "NTc.Cm2", "NTc.Cm4", "cxb1.1")
  
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
  rm(ck, ck.cc, ck.cxc, p.cc, p.cxc)
  
  # (1.2) RECEPTOR -------------------------------------------------------------
  
  # load tier2 sequence positions
  ckr.pos <- c("class", "gn1x23", "gn1x24", "gn6x58", "gn6x62", "gn7x27")
  
  # load alignment df
  ckr <- read_csv("03_ckr_seq/data/processed/ALL_classa_df.csv")
  ckr <- ckr[(names(ckr) %in% ckr.pos)]
  rm(ckr.pos)
  
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
  rm(ckr, ckr.cc, ckr.cxc, p.cc, p.cxc)
  