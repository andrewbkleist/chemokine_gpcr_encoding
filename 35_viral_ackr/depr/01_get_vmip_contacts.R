# Name:     01_get_vmip_contacts.R
# Updated:  20230106
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: FIND 4RWS (VMIP-CXCR4) OVERLAPPING CONTACTS ###########################

  # import, clean, unite
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    #filter(no_pdb >= 2) %>%
    select(file, source_gnccn, target_gnccn) 
  
  # import CC/CXC
  cc.cxc <- read_csv("32_rin_tier2/output/cc_cxc_specific_all.csv") %>% 
    filter(type != "specific") %>% select(source_gnccn, target_gnccn)
  
  # similarities and differences among complexes
  pdb.5uiw <- rin %>% filter(file == "5uiw") %>% select(-file)
  pdb.6wwz <- rin %>% filter(file == "6wwz") %>% select(-file)
  pdb.4rws <- rin %>% filter(file == "4rws") %>% select(-file)
  pdb.6lfo <- rin %>% filter(file == "6lfo") %>% select(-file)
  pdb.ngo <- rin %>% filter(file == "ngo") %>% select(-file)
  pdb.4xt1 <- rin %>% filter(file == "4xt1") %>% select(-file)
  pdb.5wb2 <- rin %>% filter(file == "5wb2") %>% select(-file)
  pdb.zheng <- rin %>% filter(file == "zheng") %>% select(-file)
  
  
  # intersection CXC - CXCL12 model, CXCL8
  intersect.cxc <- intersect(pdb.6lfo, pdb.ngo) # 15 (common to all CXC)
  
  # intersection CC - CCL20, CCL5, Zheng model
  intersect.cc <- intersect(pdb.5uiw, pdb.6wwz)
  intersect.cc <- intersect(intersect.cc, pdb.zheng) # 13 (common to all)
  
  # identify CC-/CXC-only contacts
  setdiff.cc <- setdiff(intersect.cc, pdb.6lfo)
  setdiff.cc <- setdiff(setdiff.cc, pdb.ngo)
  
  setdiff.cxc <- setdiff(intersect.cxc, pdb.5uiw)
  setdiff.cxc <- setdiff(setdiff.cxc, pdb.6wwz)
  setdiff.cxc <- setdiff(setdiff.cxc, pdb.zheng)
  
  # intsersect
  intersect.cc.cxc <- intersect(intersect.cc, intersect.cxc)
  
  # vMIP
  intersect.cc.vmip  <- intersect(setdiff.cc, pdb.4rws) # 8 (more CC-like)
  intersect.cxc.vmip  <- intersect(setdiff.cxc, pdb.4rws) # 4
  
  intersect.cc.cxc.vmip <- intersect(intersect.cc.cxc, pdb.4rws)
  
  # remove
  rm(pdb.5uiw, pdb.6wwz, pdb.4rws, pdb.6lfo, pdb.ngo, pdb.4xt1, pdb.5wb2, pdb.zheng)
