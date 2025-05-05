# Name:     01_cx3cl1_comp.R
# Updated:  20201201
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### COMPARISON OF CX3CL1 N-TERMINAL CONTACTS W US28 ##########################

  # import contacts
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(ck == "cx3cl1") %>% filter(dom1 == "NTc") %>%
    select(file, source_gnccn, target_gnccn)
  
  # subset
  pdb.4xt1 <- rin %>% filter(file == "4xt1") %>% select(-file)
  pdb.5wb2 <- rin %>% filter(file == "5wb2") %>% select(-file)
  
  # similarities, differences
  only.4xt1 <- setdiff(pdb.4xt1, pdb.5wb2) # 22!
  only.5wb2 <- setdiff(pdb.5wb2, pdb.4xt1) # 26!
  overlap <- intersect(pdb.4xt1, pdb.5wb2)
  
  # add labels, combine
  only.4xt1$type <- c("4xt1_only")
  only.5wb2$type <- c("5wb2_only")
  overlap$type <- c("overlap")
  master <- bind_rows(only.4xt1, only.5wb2, overlap)
  
  # remove
  rm(pdb.4xt1, pdb.5wb2, only.4xt1, only.5wb2, overlap)
  
  
  # NTc agnostic (receptor contacts only)
  pdb.4xt1.ag <- rin %>% filter(file == "4xt1") %>% select(-file, -source_gnccn) %>% unique()
  pdb.5wb2.ag <- rin %>% filter(file == "5wb2") %>% select(-file, -source_gnccn) %>% unique()
  
  only.4xt1.ag <- setdiff(pdb.4xt1.ag, pdb.5wb2.ag) # 9!
  only.5wb2.ag <- setdiff(pdb.5wb2.ag, pdb.4xt1.ag) # 11!
  overlap.ag <- intersect(pdb.4xt1.ag, pdb.5wb2.ag) # 9!
  
  # add labels, combine
  only.4xt1.ag$type <- c("4xt1_only")
  only.5wb2.ag$type <- c("5wb2_only")
  overlap.ag$type <- c("overlap")
  master.ag <- bind_rows(only.4xt1.ag, only.5wb2.ag, overlap.ag)
  
  # remove
  rm(pdb.4xt1.ag, pdb.5wb2.ag, only.4xt1.ag, only.5wb2.ag, overlap.ag)
 