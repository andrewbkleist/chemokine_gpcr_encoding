# Name:     01_motif_heatmaps.R
# Updated:  20210518
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### 1: HEATMAP - CHEMOKINE ###################################################
   
  # (1.0) HOW MANY SLiMs? --------------------------------------------------------
  # count how many unique fragments there are
  data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")  %>%
    select(motif) %>%
    unique()
  rm(data)
  
  # count how many fragments are conserved among orthologs
  data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")  %>%
    filter(pct_ortho >= 0.5) %>%
    select(motif) %>%
    unique()
  rm(data)
  
  # count how many fragments are are found in > 5 receptors
  data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")  %>%
    filter(count_super >= 5) %>%
    filter(pct_ortho >= 0.5) %>%
    select(motif) %>%
    unique()
  rm(data)
  
  # count how many fragments are are found in 1 receptors
  data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")  %>%
    filter(count_super == 1) %>%
    filter(pct_ortho >= 0.5) %>%
    select(motif) %>%
    unique()
  rm(data)

  # (1.1) IMPORT, ORDER --------------------------------------------------------

    # import, reformat
    data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")  
  
    # select only relevant info (motif, protein, class, pct ortho)
    data <- data %>% select(motif, protein, class, pct_ortho)
    
    # order chemokines
    order.ck <- as.factor(tolower(c("ccl1","ccl2", "ccl3", "ccl3l1",
                                    "ccl4","ccl4l1","ccl5","ccl7","ccl8",
                                    "ccl11","ccl13","ccl14","ccl15",
                                    "ccl16","ccl17","ccl18","ccl19",
                                    "ccl20","ccl21","ccl22","ccl23",
                                    "ccl24","ccl25","ccl26","ccl27",
                                    "ccl28","cxcl1","cxcl2","cxcl3",
                                    "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                    "cxcl7","cxcl8","cxcl9","cxcl10",
                                    "cxcl11","cxcl12","cxcl13", "cxcl14", "cxcl16",
                                    "cxcl17", "cx3cl1", "xcl1", "xcl2")))
    levels(data$protein)
    data$protein <- factor(data$protein, levels = order.ck)
    
    # order motifs based on appearance:
    # want most unique to most conserved (global) then most to least conserved among orthologs
    data <- data %>% add_count(motif)
    data <- data %>% arrange(n, protein, -pct_ortho)
    order.motif <- unique(data$motif)
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.motif))
    
  #(1.2) PLOTTING --------------------------------------------------------------
    
    # MASTER HEATMAP
    data %>%
      # filter(motif == "AA") %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4")+
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    
##### 2: HEATMAP - RECEPTOR NT #################################################
    
    # (2.0) HOW MANY SLiMs? --------------------------------------------------------
    # count how many unique motifs there are
    data <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")  %>%
      select(motif) %>%
      unique()
    rm(data)
    
    # count how many motifs are conserved among orthologs
    data <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")  %>%
      filter(pct_ortho >= 0.5) %>%
      select(motif) %>%
      unique()
    rm(data)
    
    # count how many motifs are are found in > 5 receptors
    data <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")  %>%
      filter(count_super >= 5) %>%
      filter(pct_ortho >= 0.5) %>%
      select(motif) %>%
      unique()
    rm(data)
    
    
    
    
    
    # (2.1) IMPORT, ORDER --------------------------------------------------------
    
    # import, reformat
    data <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")  
    
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
    
    # order motifs based on appearance:
    # want most unique to most conserved (global) then most to least conserved among orthologs
    data <- data %>% add_count(motif)
    data <- data %>% arrange(n, protein, -pct_ortho)
    order.motif <- unique(data$motif)
    levels(data$motif)
    data$motif <- factor(data$motif, levels = rev(order.motif))
    
    
    #(2.2) PLOTTING ----------------------------------------------------------------
    
    # MASTER HEATMAP
    # use below to find position of motif; must move "fill" inside aes to work
    # data <- data %>% mutate(dyd = case_when(
    #   motif == "DYD" ~ "yes",
    #   motif != "DYD" ~ "no"
    # ))
    # data$dyd <- as.factor(data$dyd)
    
    data %>%
      ggplot() + 
      geom_tile(aes(protein, motif), fill = "mediumorchid4") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

 