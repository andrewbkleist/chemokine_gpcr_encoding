# Name:     ***
# Updated:  ***
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
rin <- rin %>% 
  filter(file != "6meo") %>%
  filter(all_para_ck <0.5 | all_non_ackr_para_ckr <0.5) %>%
  # filter(source_gnccn == "cxb1.1" & target_gnccn == "7x27") %>%
  count(source_gnccn, target_gnccn) %>%
  filter(n > 7)
  



rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
rin <- rin %>% 
  filter(file != "6meo") %>%
  filter(ortho_cons_ck >=0.5 & ortho_cons_ckr >=0.5) %>%
  filter(all_para_ck <0.5 & all_para_ckr <0.5) %>%
  select(file, ck, ckr, source_gnccn, target_gnccn, 
         all_para_ck, all_non_ackr_para_ckr, 
         cc_cxc_lr_ck, cc_cxc_lr_ckr,
         ortho_cons_ck, ortho_cons_ckr) %>%
  count(source_gnccn, target_gnccn)
  
  