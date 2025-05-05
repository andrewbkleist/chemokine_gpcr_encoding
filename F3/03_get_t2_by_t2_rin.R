# Name:     02_get_tier2_rin.R
# Updated:  20221111
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE (20221111): CC & CXC COMPLEXES
# (1) CCL5 : CCR5 & VARIANTS (SEMI-DEGENERATE)
# 5uiw: ccl5-5p7 / ccr5
# 7o7f: ccl5-6p4 / ccr5 **
# 7f1r: ccl5 / ccr5 (incomplete) **
# zheng: ccl5 / ccr5 (model)
# 
# (2) OTHER CCL* : CCR* COMPLEXES
# 7vl9: ccl15 / ccr1
# 7xa3: ccl2 / ccr2
# 7f1t: ccl3 / ccr5 
# 6wwz: ccl20 / ccr6
#
# (3) CXCL* : CXCR* COMPLEXES
# 6lfo: cxcl8 / cxcr2
# ngo: cxcl12 / cxcr4  (model)
#
# (4) EXCLUDES...
# 7sk3: cxcl12 / ackr3
# 7xbx: cx3cl1 / cx3cr1
# 4rws: vmipii / cxcr4
# 4xt1: cx3cl1 / us28
# 5wb2: cx3cl1.35 / us28 **
# 6meo: gp120 / ccr5 (HIV complex)

################################################################################
  
  # import data, select only CC/CXC
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                       "7xa3", "7f1t", "6wwz", "6lfo", "ngo"))
  
  
  # test <- rin %>% filter(cc_cxc_lr_ck >= 0.75, cc_cxc_lr_ckr >= 0.75) %>%
  #   select(source_gnccn, target_gnccn) %>% unique()
  
  # select tier 2 x tier 2 contacts (175)
  rin <- rin %>% filter(cc_cxc_lr_ck >= 0.75, cc_cxc_lr_ckr >= 0.75) %>%
    filter(! source_gnccn %in% c("NTc.Cm10", "NTc.Cm9", "NTc.Cm8","NTc.Cm7","NTc.Cm6",
                             "NTc.Cm5"))
  
  # add cc, cxc information
  rin <- rin %>% mutate(cc_cxc = case_when(
    grepl("cxcl",ck) ~ "cxc",
    grepl("ccl",ck) ~ "cc"
  ))
  
  # RESIDUE LEVEL - count number of unique contacts
  # first need to make CCL5-CCR5 degenerate, then remove redundancies
  rin <- rin %>% mutate(file = case_when(
    ck == "ccl5" ~ "ccl5_ccr5",
    ck != "ccl5" ~ file
  ))
  rin <- rin %>% select(file:dom2, cc_cxc) %>% unique()
  rin.unique <- rin %>% count(source_gnccn, target_gnccn, cc_cxc) %>% unique()
  rin.unique <- rin.unique %>% pivot_wider(names_from = cc_cxc, values_from = n)
  rin.unique[is.na(rin.unique)] <- 0
  
  # RESIDUE LEVEL - identify CC, CXC consensus
  cc.cons <- rin.unique %>% filter(cc > 2 & cxc == 0)
  cxc.cons <- rin.unique %>% filter(cc ==0 & cxc == 2)
  shared.cons <- rin.unique %>% filter(cc >2 & cxc == 2)
  
    # RESULT:
    # (1) CC-specific consensus contacts: 
    # - NTc.Cm3	1x24
    # - NTc.Cm3	1x28
    # - NTc.Cm4	1x28
    # - b1b2.6	45x51
    # - b1b2.9	5x36
    # - b1b2.10	5x36
    # - (also b1b2.6 ECL2.Cp6)
    # (2) CXC-specific consensus contacts: 
    # - NTc.Cm1	6x62
    # - NTc.Cm1	7x27
    # - NTc.Cm3	6x58
    # - NTc.Cm3	7x34
    # (3) Shared consensus contacts:
    # - cxb1.1  7x27
  
  # RESIDUE LEVEL - isolate ALL unique ligand, receptor residues among CC, CXC
  ck.count <- rin %>% count(source_gnccn, cc_cxc)
  ck.count <- ck.count %>% pivot_wider(names_from = cc_cxc, values_from = n)
  ck.count[is.na(ck.count)] <- 0
  
  ckr.count <- rin %>% count(target_gnccn, cc_cxc)
  ckr.count <- ckr.count %>% pivot_wider(names_from = cc_cxc, values_from = n)
  ckr.count[is.na(ckr.count)] <- 0
  
    # RESULT:
    # (1) CC chemokine hotspots: Ntc.Cm3, NTc.Cm4
    # (2) CXC chemokine hotspots: Ntc.Cm1, NTc.Cm3
    # (3) CC receptor hotspots: 1x28, 1x24, 5x36, 7x31
    # (4) CXC receptor hotspots: 45x51, 7x27, 7x34, 6x58, 6x62
  
  # DOMAIN LEVEL - count unique domain arrangements among CC, CXC PAIRWISE
  rin.dom <- rin %>% count(dom1, dom2, cc_cxc) %>% unique()
  rin.dom <- rin.dom %>% pivot_wider(names_from = cc_cxc, values_from = n)
  rin.dom[is.na(rin.dom)] <- 0
  
    # RESULT: 
    # (1) CCL-CCR complexes use NTc-TM1, NTc-TM7, b1b2-ECL2, b1b2-TM5 
    # (2) CXCL-CXCR complexes use NTc-TM6, NTc-TM7, NTc-ECL2
  
  # DOMAIN LEVEL - count unique domain arrangements among CC, CXC - INDIVIDUAL
  rin.dom.ck <- rin %>% count(dom1, cc_cxc)
  rin.dom.ck <- rin.dom.ck %>% pivot_wider(names_from = cc_cxc, values_from = n)
  rin.dom.ck[is.na(rin.dom.ck)] <- 0
  
  rin.dom.ckr <- rin %>% count(dom2, cc_cxc)
  rin.dom.ckr <- rin.dom.ckr %>% pivot_wider(names_from = cc_cxc, values_from = n)
  rin.dom.ckr[is.na(rin.dom.ckr)] <- 0
  
    # RESULT: 
    # (1) Subfamily specific residues for CC on chemokine are Ntc then b1b2
    # (2) Subfamily specific residues for CXC on chemokine are Ntc
    # (3) Subfamily specific residues for CC on receptor are TM1, TM7
    # (4) Subfamily specific residues for CXC on receptor are TM6, TM7, ECL2
  
  
##### NO FILTERING BASED ON T2 - CONSERVED CONTACTS ############################
  # import data, select only CC/CXC
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                       "7xa3", "7f1t", "6wwz", "6lfo", "ngo"))
  
  # select tier 2 x tier 2 contacts (175)
  # rin <- rin %>% filter(cc_cxc_lr_ck >= 0.75, cc_cxc_lr_ckr >= 0.75) %>%
  #   filter(! source_gnccn %in% c("NTc.Cm10", "NTc.Cm9", "NTc.Cm8","NTc.Cm7","NTc.Cm6",
  #                                "NTc.Cm5"))
  
  # add cc, cxc information
  rin <- rin %>% mutate(cc_cxc = case_when(
    grepl("cxcl",ck) ~ "cxc",
    grepl("ccl",ck) ~ "cc"
  ))
  
  # RESIDUE LEVEL - count number of unique contacts
  # first need to make CCL5-CCR5 degenerate, then remove redundancies
  rin <- rin %>% mutate(file = case_when(
    ck == "ccl5" ~ "ccl5_ccr5",
    ck != "ccl5" ~ file
  ))
  rin <- rin %>% select(file:dom2, cc_cxc) %>% unique()
  
  rin.unique <- rin %>% count(source_gnccn, target_gnccn, cc_cxc) %>% unique()
  rin.unique <- rin.unique %>% pivot_wider(names_from = cc_cxc, values_from = n)
  rin.unique[is.na(rin.unique)] <- 0
  
  
##### T1-T2 CONTACTS #################################################
  

  # import data, select only CC/CXC
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                       "7xa3", "7f1t", "6wwz", "6lfo", "ngo"))
  
  # select tier 2 x tier 2 contacts (175)
  rin1 <- rin %>% filter(all_para_ck >= 0.75, cc_cxc_lr_ckr >= 0.75) %>%
    filter(! source_gnccn %in% c("NTc.Cm10", "NTc.Cm9", "NTc.Cm8","NTc.Cm7","NTc.Cm6",
                                 "NTc.Cm5"))
  
  rin2 <- rin %>% filter(cc_cxc_lr_ck >= 0.75, all_non_ackr_para_ckr >= 0.75) %>%
    filter(! source_gnccn %in% c("NTc.Cm10", "NTc.Cm9", "NTc.Cm8","NTc.Cm7","NTc.Cm6",
                                 "NTc.Cm5"))
  
  rin <- rbind(rin1, rin2)
  
  # add cc, cxc information
  rin <- rin %>% mutate(cc_cxc = case_when(
    grepl("cxcl",ck) ~ "cxc",
    grepl("ccl",ck) ~ "cc"
  ))
  
  # RESIDUE LEVEL - count number of unique contacts
  # first need to make CCL5-CCR5 degenerate, then remove redundancies
  rin <- rin %>% mutate(file = case_when(
    ck == "ccl5" ~ "ccl5_ccr5",
    ck != "ccl5" ~ file
  ))
  rin <- rin %>% select(file:dom2, cc_cxc) %>% unique()
  
  rin.unique <- rin %>% count(source_gnccn, target_gnccn, cc_cxc) %>% unique()
  rin.unique <- rin.unique %>% pivot_wider(names_from = cc_cxc, values_from = n)
  rin.unique[is.na(rin.unique)] <- 0
  
  # RESIDUE LEVEL - identify CC, CXC consensus
  cc.cons <- rin.unique %>% filter(cc > 2 & cxc == 0)
  cxc.cons <- rin.unique %>% filter(cc ==0 & cxc == 2)
  shared.cons <- rin.unique %>% filter(cc >2 & cxc == 2)
  