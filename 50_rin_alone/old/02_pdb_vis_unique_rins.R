# Name:     02_pdb_vis_unique_contacts.R
# Updated:  20210408
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: GET RESNOS FOR UNIQUE AND CONSERVED STRUCTURAL CONTACTS ###############
  
  # import
  rin.proc <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  rin.raw <- read_csv("01_structure_contacts/output/RIN_residue.csv") %>%
    filter(class == "full") %>%
    select(ResNum1, ResNum2, source_gnccn, target_gnccn, file) %>% unique()
  rin.proc <- left_join(rin.proc, rin.raw)
  rm(rin.raw)
  
  rin.proc <- rin.proc %>% select(file:target_gnccn, ResNum1, ResNum2, dom1:cancer_mut_count_ckr)
  
  
  