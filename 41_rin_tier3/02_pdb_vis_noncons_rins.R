# Name:     02_pdb_vis_noncons_contacts.R
# Updated:  20210414
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: GET RESNOS FOR NONCONSERVED STRUCTURAL CONTACTS #######################
  
  # import
  rin.proc <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  rin.raw <- read_csv("01_structure_contacts/output/RIN_residue.csv") %>%
    filter(class == "full") %>%
    select(ResNum1, ResNum2, source_gnccn, target_gnccn, file) %>% unique()
  rin.proc <- left_join(rin.proc, rin.raw)
  rm(rin.raw)
  
  # reorder
  rin.proc <- rin.proc %>%
    select(file:target_gnccn, ResNum1, ResNum2, dom1:cancer_mut_count_ckr)
  
  # subset
  rin.proc.cons <- rin.proc %>% filter( (all_para_ck < 0.5) & (all_non_ackr_para_ckr < 0.5) ) %>% # PARALOG: 348 total noncons-noncons
    filter((ortho_cons_ck > 0.5) & (ortho_cons_ckr > 0.5))#  %>% # ORTHOLOG: 112 cons-cons, 12 noncons-noncons, 93 cons-non; 4 non-cons; remainder NA
  #filter(cc_cxc_lr_ck < 0.75 | cc_cxc_lr_ckr < 0.75)

  
  
  
  