source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
# (1) SUBFAMILY-PREDICTIVE BY SUBFAMILY-PREDICTIVE ONLY ------------------------
# import data, select only CC/CXC
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0"))
  
# select subfamily x subfamily contacts
rin <- rin %>% filter(cc_cxc_lr_ck >= 0.75, cc_cxc_lr_ckr >= 0.75)

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
rin.unique <- rin %>% dplyr::count(source_gnccn, target_gnccn, cc_cxc) %>% unique()
rin.unique <- rin.unique %>% pivot_wider(names_from = cc_cxc, values_from = n)
rin.unique[is.na(rin.unique)] <- 0
  
# RESIDUE LEVEL - identify CC, CXC consensus
cc.cons <- rin.unique %>% filter(cc > 2 & cxc == 0)
cxc.cons <- rin.unique %>% filter(cc ==0 &  cxc > 1)
shared.cons <- rin.unique %>% filter(cc >2 &  cxc > 1)

  # RESULT:
  # (1) CC-specific consensus contacts: 
  cc.cons
  # - B1.7         ECL2.Cp3
  # - NTc.Cm3      1x24
  # - NTc.Cm3      1x28
  # - NTc.Cm4      1x28
  # - b1b2.10      5x36
  # - b1b2.16      5x32
  # - b1b2.16      ECL2.Cp6
  # - b1b2.4       ECL2.Cp3
  # - b1b2.6       45x51
  # - b1b2.9       5x32
  # - b1b2.9       5x36
  # (2) CXC-specific consensus contacts: 
  cxc.cons
  # - NTc.Cm1      6x58
  # - NTc.Cm1      7x27
  # - NTc.Cm3      6x58
  # - NTc.Cm3      7x34
  # (3) Shared consensus contacts:
  shared.cons  
  # cxb1.1 7x27

# (2) SUBFAMILY-PREDICTIVE BY CONSERVED ----------------------------------------
# import data, select only CC/CXC
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0"))

# select tier 2 x tier 2 contacts (175)
rin1 <- rin %>% filter(all_para_ck >= 0.75, cc_cxc_lr_ckr >= 0.75)
rin2 <- rin %>% filter(cc_cxc_lr_ck >= 0.75, all_non_ackr_para_ckr >= 0.75) 
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
rin <- rin %>% dplyr::select(file:dom2, cc_cxc) %>% unique()

rin.unique <- rin %>% dplyr::count(source_gnccn, target_gnccn, cc_cxc) %>% unique()
rin.unique <- rin.unique %>% pivot_wider(names_from = cc_cxc, values_from = n)
rin.unique[is.na(rin.unique)] <- 0
  
# RESIDUE LEVEL - identify CC, CXC consensus
cc.cons <- rin.unique %>% filter(cc > 2 & cxc == 0)
cxc.cons <- rin.unique %>% filter(cc ==0 & cxc > 1)
shared.cons <- rin.unique %>% filter(cc >2 & cxc > 1)

  # RESULT:
  # (1) CC-specific consensus contacts: 
  cc.cons
  # - NTc.Cm1      1x22
  # - b1b2.12      5x32
  # - b1b2.12      5x36
  # - b3h.2        ECL2.Cp6
  # (2) CXC-specific consensus contacts: 
  cxc.cons
  # - CX.1 7x27
  # - CX.4 1x22
  # (3) Shared consensus contacts:
  shared.cons  
  # (none)
  