# Name:     01_ck_ckr_ridge_plots.R
# Updated:  20230106
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggridges)

##### 1: CHEMOKINE RIDGE PLOT ##################################################

  # (1.1) Get CK (i) T2 positions (ii) @ interface -----------------------------

    # import T2
    ck.t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      filter(mean >= 0.75)
    
    # import interface
    ck.rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
      select(source_gnccn) %>% unique()
    ck.rin <- ck.rin$source_gnccn
    
    # subset T2 for those at interface...
    ck.t2 <- ck.t2 %>% filter(motif %in% ck.rin)
    ck.t2 <- ck.t2 %>% filter(motif != "NTc.Cm10" & motif != "NTc.Cm9" & motif != "NTc.Cm8" & motif != "NTc.Cm7")
    ck.t2 <- ck.t2$motif
    
    # import chemokine scores, subset by T2 at interface
    ck.paralog.acc <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
    ck.virus.acc <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")    
    
    ck.paralog.acc <- ck.paralog.acc %>% filter(position %in% ck.t2)
    ck.virus.acc <- ck.virus.acc %>% filter(position %in% ck.t2)
    
    # remove
    rm(ck.t2, ckr.t2)
    rm(ck.rin, ckr.rin)
  
  # (1.2) CHEMOKINE - SUBSET PROTEINS OF INTEREST AND GRAPH --------------------
  
    # subset, combine
    viral.sub <- c("k6xhhv8p", "vmip2xhhv8p", "vcxcl1xhcmvm") #could get vCXCL4a/b, add vccl3xhhv8p
    ck.virus.acc <- ck.virus.acc %>% filter(protein %in% viral.sub)
    
    ck.paralog.acc <- bind_rows(ck.paralog.acc, ck.virus.acc)
    rm(viral.sub, ck.virus.acc)
    
    # order
    order <- c("ccl1","ccl2","ccl3","ccl3l1","ccl4",
               "ccl4l1","ccl5","ccl7","ccl8","ccl11",
               "ccl13","ccl14","ccl15","ccl16","ccl17",
               "ccl18","ccl19","ccl20","ccl21","ccl22",
               "ccl23","ccl24","ccl25","ccl26","ccl27",
               "ccl28","cxcl1","cxcl2","cxcl3","cxcl4",
               "cxcl4l1","cxcl5","cxcl6","cxcl7","cxcl8",
               "cxcl9","cxcl10","cxcl11","cxcl12","cxcl13",
               "cxcl14","cxcl16","cxcl17","cx3cl1","xcl1",
               "xcl2","k6xhhv8p","vmip2xhhv8p","vcxcl1xhcmvm")
    
    ck.paralog.acc$protein <- factor(ck.paralog.acc$protein, levels = (order))
    
    # plot all
    ck.paralog.acc %>%
      select(-score) %>%
      unique() %>%
      ggplot(aes(x = mean, y = protein)) +
      geom_density_ridges(fill = "steelblue3") +
      # geom_hline(yintercept = 0.6) +
      # geom_hline(yintercept = 0.4) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0)) +
      coord_flip()
    
    # plot viral only
    ck.paralog.acc %>%
      select(-score) %>%
      filter(protein == "vmip2xhhv8p" | protein == "k6xhhv8p" | protein == "vcxcl1xhcmvm") %>%
      unique() %>%
      ggplot(aes(x = mean, y = protein)) +
      geom_density_ridges(fill = "steelblue3") +
      theme_minimal()
  
    # remove
    rm(ck.paralog.acc, order)

##### 2: RECEPTOR RIDGE PLOT ###################################################
    
  # (2.1) Get CKR (i) T2 positions (ii) @ interface ----------------------------
    
    # import T2
    ckr.t2 <- read_csv("03_ckr_seq/output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      filter(mean > 0.75)
    
    # import interface
    ckr.rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
      select(target_gnccn) %>% unique()
    ckr.rin$gn <- c("gn")
    ckr.rin <- ckr.rin %>% unite(col = target_gnccn, c(gn, target_gnccn), sep = "", remove = TRUE)
    ckr.rin <- ckr.rin$target_gnccn
    
    # subset T2 for those at interface...
    ckr.t2 <- ckr.t2 %>% filter(motif %in% ckr.rin)
    ckr.t2 <- ckr.t2 %>% filter(motif != "gnNTr.Cm26" & motif != "gnNTr.Cm11" & motif != "gnNTr.Cm10" & motif != "gnNTr.Cm6" & motif != "gnNTr.CmCp6")
    ckr.t2 <- ckr.t2$motif
    
    # import receptor scores
    ckr.classa.acc <- read_csv("03_ckr_seq/output/CKR_CLASSA_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
    colnames(ckr.classa.acc)[1] <- c("score")
    ckr.virus.acc <- read_csv("03_ckr_seq/output/CKR_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
    
    ckr.classa.acc <- ckr.classa.acc %>% filter(position %in% ckr.t2)
    ckr.virus.acc <- ckr.virus.acc %>% filter(position %in% ckr.t2)
  
    # remove
    rm(ckr.t2)
    rm(ckr.rin)
    
  # (2.2) RECEPTOR - SUBSET PROTEINS OF INTEREST AND GRAPH ---------------------
  
    # subset, combine - CLASSA
    classa.sub <- c("ccr1","ccr2", "ccr3", "ccr4",
                    "ccr5","ccr6","ccr7","ccr8","ccr9",
                    "ccr10","cxcr1","cxcr2","cxcr3",
                    "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                    "ackr1","ackr2","ackr3","ackr4",
                    "ccrl2")
    ckr.classa.acc <- ckr.classa.acc %>% filter(protein %in% classa.sub)
    
    
    # subset, combine - VIRAL
    viral.sub <- c("us28xhcmva", "vu51xhhv6u") # missing U12 HHV7 and HHV6 (different sets of CKs)
    ckr.virus.acc <- ckr.virus.acc %>% filter(protein %in% viral.sub)
    
    ckr.classa.acc <- bind_rows(ckr.classa.acc, ckr.virus.acc)
    rm(viral.sub, ckr.virus.acc)
  
    # order
    order <- c("ccr1","ccr2", "ccr3", "ccr4",
               "ccr5","ccr6","ccr7","ccr8","ccr9",
               "ccr10","cxcr1","cxcr2","cxcr3",
               "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
               "ackr1","ackr2","ackr3","ackr4",
               "ccrl2", "us28xhcmva", "vu51xhhv6u")
    
    ckr.classa.acc$protein <- factor(ckr.classa.acc$protein, levels = (order))
    
    # plot
    ckr.classa.acc %>%
      # select(-score) %>%
      unique() %>%
      ggplot(aes(x = mean, y = protein)) +
      geom_density_ridges(fill = "steelblue3") +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0)) +
      coord_flip()
    
    ckr.classa.acc %>%
      filter(protein == "ackr1" | protein == "ackr2" | protein == "ackr3" | protein == "ackr4" | protein == "ccrl2" | protein == "us28xhcmva" | protein == "vu51xhhv6u") %>%
      unique() %>%
      ggplot(aes(x = mean, y = protein)) +
      geom_density_ridges(fill = "steelblue3") +
      theme_minimal() 
    
    # remove
    rm(classa.sub, ckr.virus.acc, order)
    
