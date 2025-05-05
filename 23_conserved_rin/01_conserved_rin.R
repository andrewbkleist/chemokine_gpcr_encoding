# Name:     01_conserved_rin.R
# Updated:  20221106
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE (20221106): 
# IDENTIFYING STRUCTURALLY CONSERVED RESIDUE-RESIDUE INTERACTIONS
# There are 15 complexes (excluding HIV-CCR5 structure), 5 degenerate 
# CCL5-CCR5 complexes (5P7, 6P7, WT, model), 4 additional CCL:CCR complexes,
# 2x CXCL:CXCR complexes, 2 noncanonical complexes (ACKR3, CX3CL1:CX3CR1),
# and 3 viral (vMIP-II:CXCR4 and 2x CX3CL1:US28).
#
# 2/3 of complexes would be >/= 10complexes 

##### 1:  GET CONSERVED  CHEMOKINE-GPCR CONTACTS ###############################

  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file != "6meo") %>%
    select(source_gnccn, target_gnccn, no_pdb, all_para_ck, all_non_ackr_para_ckr) %>% 
    unique() %>%
    filter(no_pdb > 7)
  
  
  
  # plot
  rin %>%
    # filter(no_pdb >4) %>%
    ggplot(aes(all_non_ackr_para_ckr, all_para_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
    # scale_size(range = c(2, 8)) +
    xlim(0,1) +
    ylim(0,1) +
    theme_minimal()
  
  
  # rin <- rin %>% filter(all_ortho_ck > 0.4 & all_ortho_ckr > 0.6)
  rin <- rin %>% filter((all_cc_cxc_para_ck > 0.5) & (all_non_ackr_para_ckr > 0.5))  %>%
     select(source_gnccn, target_gnccn, no_pdb) %>% unique()
  
  
  