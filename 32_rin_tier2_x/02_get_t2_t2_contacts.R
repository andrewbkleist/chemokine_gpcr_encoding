# Name:     02_get_t2_t2_contacts.R
# Updated:  20201105
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)

##### 1: IDENTIFY TIER 2 TO TIER 2 CONTACTS ####################################

  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  rin.con <- read_csv("20_rin_tier1/output/rin_consensus.csv")
  rin.cc.cxc <- read_csv("32_rin_tier2/output/cc_cxc_specific_cyto.csv")
  rin.sum <- bind_rows(rin.con, rin.cc.cxc)
  rm(rin.con, rin.cc.cxc)
  rin <- left_join(rin, rin.sum)  
  #rm(rin.sum)
  
  # get contacts > 0.75 AND "key contact" (consensus, cc or cxc) or others that are found...
  # exclude CX3CL1-US28 complexes
  rin.a <- rin %>% filter(cc_cxc_lr_ck >= 0.8 & cc_cxc_lr_ckr >= 0.8) %>% filter(ckr != "us28")
  rin.b <- rin %>% filter(!is.na(type)) %>%  filter(cc_cxc_lr_ck > 0.75 & cc_cxc_lr_ckr > 0.75)
  rin <- bind_rows(rin.a, rin.b)
  rm(rin.a, rin.b)
  rin.t2 <- rin %>%
    select(source_gnccn, target_gnccn, no_pdb) %>%
    unique()
  
  # clean up; remove distal N-terminal 
  rin.t2 <- left_join(rin.t2, rin.sum)
  # rin.t2 <- rin.t2 %>% filter(source_gnccn != "NTc.Cm8") %>% filter(source_gnccn != "NTc.Cm9") %>%
  #   filter(target_gnccn != "NTr.Cm10") %>% filter(target_gnccn != "ECL2.Cp6")
  rin.t2.ck <- rin.t2 %>% select(source_gnccn) %>% unique()
  rin.t2.ckr <- rin.t2 %>% select(target_gnccn) %>% unique()
  rm(rin.sum)
  
  rin.t2 <- rin.t2 %>% mutate(type=case_when(
    is.na(type) ~ "specific",
    !is.na(type) ~ type
  ))
  
  # write
  write_csv(rin.t2, "32_rin_tier2/output/cc_cxc_specific_all.csv")
  
  # remove
  rm(rin, rin.t2, rin.t2.ck, rin.t2.ckr)
  
