# Name:     01_get_t2_logo.R
# Updated:  20201113
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

##### 1: TIER 2 SEQUENCE LOGOS #################################################

  # (1.1) CHEMOKINE- CONSENSUS -------------------------------------------------
    
    # load tier2 sequence positions
    ck.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_all.csv") %>% 
      filter(type == "consensus") %>%
      select(source_gnccn)
    
    ck.pos <- c("class", ck.pos$source_gnccn)
    
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
  
  # (1.2) RECEPTOR - CONSENSUS -------------------------------------------------------------
  
    # load tier2 sequence positions
    ckr.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_all.csv") %>%
      filter(type == "consensus") %>%
      select(target_gnccn)
    ckr.pos$gn <- c("gn")
    ckr.pos <- ckr.pos %>% unite(target_gnccn, c(gn, target_gnccn), sep = "")
    ckr.pos <- c("class", ckr.pos$target_gnccn)
    
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
    
    
##### TEMP TEMP TEMP 1: TIER 2 SEQUENCE LOGOS #################################################
    
    # (1.1) CHEMOKINE ------------------------------------------------------------
    
    # load CC and CXC-specific contacts
    ck.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_cyto.csv") %>% select(source_gnccn)
    ck.pos <- c("class", ck.pos$source_gnccn)
    
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
    
    # load CC and CXC-specific contacts
    ckr.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_cyto.csv") %>% select(target_gnccn)
    ckr.pos$gn <- c("gn")
    ckr.pos <- ckr.pos %>% unite(target_gnccn, c(gn, target_gnccn), sep = "")
    ckr.pos <- c("class", ckr.pos$target_gnccn, "gnECL3.5") 
    
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
    
    
    # FUNCTION IN PROGRESS
    
    GetLogo <- function(RINTYPE, SOURCE.TARGET.GNCCN, ALN){
      
      # load tier2 sequence positions
      ck.pos <- read_csv("32_rin_tier2/output/cc_cxc_specific_all.csv") %>% 
        filter(type == RINTYPE) %>%
        select(source_gnccn)
      
      #ck.pos <- c("class", ck.pos$source_gnccn)
      cmd <- noquote(paste0("ck.pos <- c("class", ck.pos$", i, ")"))
      eval(parse(text = cmd))
      
      # load alignment df
      ck <- read_csv(ALN)
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
    }
    
  
  