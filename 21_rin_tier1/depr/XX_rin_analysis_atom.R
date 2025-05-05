# Name:     10_rin_analysis.R
# Updated:  20210323
# Author:   Andrew Kleist
# Figures:  Makes (1) master contact matrix and 
#           (2) histogram of no. PDBs per contact 

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
  
##### 1: GET ALL ATOMIC CONTACTS FOR T1 RESIDUES - 0.5 A CUTOFF  ###############
  
  # import rin
  rin <- read_csv("01_structure_contacts/output/RIN_atom.csv") %>% filter(class == "full") %>% filter(Chain1 != Chain2)
  # rin <- read_csv("01_structure_contacts/output/RIN_atom_1angstrom.csv") %>% filter(class == "full") %>% filter(Chain1 != Chain2)
  rin <- rin %>%
    select(source_gnccn, target_gnccn, Chain.Types, Number.of.atomic.contacts, file)
  chain.type.cnt <- rin %>% count(source_gnccn, target_gnccn, Chain.Types, file)
  
  # select B3.3
  rin.B3.3 <- rin %>% 
    filter(source_gnccn == "cxb1.2") %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts,Chain.Types, file) %>% 
    unique()
  
  rin.B3.3 <- rin.B3.3 %>%
    mutate(partner_contact = case_when(
      Chain.Types == "S-M" ~ "mc",
      Chain.Types == "M-M" ~ "mc",
      Chain.Types == "S-S" ~ "sc",
      Chain.Types == "M-S" ~ "sc"
    ))
  
  rin.B3.3 <- rin.B3.3 %>% group_by(partner_contact, file) %>%
    mutate(sum = sum(Number.of.atomic.contacts)) %>% ungroup()
    
  rin.B3.3 <- rin.B3.3 %>% select(source_gnccn, target_gnccn, partner_contact, file, sum) %>% unique()
  
  rin.B3.3 %>%
    ggplot(aes(partner_contact, sum)) +
    geom_boxplot() +
    #geom_violin(trim = F) +
    theme_minimal()
  
  stat.test <- t.test(sum ~ partner_contact, data = rin.B3.3, exact = FALSE)
  stat.test$p.value
  
  
  rin.B3.3 <- rin.B3.3 %>% group_by(partner_contact) %>% summarise(mean = mean(sum), sd = sd(sum)) %>% ungroup()
  
  rin.B3.3 %>%
    ggplot(aes(partner_contact, mean)) +
    geom_bar(stat = "identity",  position=position_dodge()) +
    geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                  position=position_dodge(.9)) + 
    theme_minimal()
  
 
  
  
  #---
  
  
  # (1) B3.3-NTr.Cm3
  rin.B3.3.NTr.Cm3 <- rin %>% 
    filter( (source_gnccn == "B3.3") & (target_gnccn == "NTr.Cm3" )) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts,Chain.Types, file) %>% 
    unique()
  rin.B3.3.NTr.Cm3 <- rin.B3.3.NTr.Cm3 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
    # 6 complexes; 6/6 are MAIN/SIDE-TO-MAIN
    # hence SIDE CHAIN or NTr.Cm3 doesn't matter
  
  # (2) B3.3-NTr.Cm1
  rin.B3.3.NTr.Cm1 <- rin %>% 
    filter( (source_gnccn == "B3.3") & (target_gnccn == "NTr.Cm1" )) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts, Chain.Types, file) %>% 
    unique()
  rin.B3.3.NTr.Cm1 <- rin.B3.3.NTr.Cm1 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
    # 6 complexes make contact; 5/6 are MAIN/SIDE-TO-MAIN; some 6/6 are MAIN/SIDE-TO-SIDE
    # hence must be compatible with side chain...
    # BUT side chain is *almost* always Pro or Leu
  
  # (3) B3.3-NTr.Cm2
  rin.B3.3.NTr.Cm2 <- rin %>% 
    filter( (source_gnccn == "B3.3") & (target_gnccn == "NTr.Cm2" )) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts, Chain.Types, file) %>% 
    unique()
  rin.B3.3.NTr.Cm2 <- rin.B3.3.NTr.Cm2 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
  # 6 complexes make contact; 5/6 are MAIN/SIDE-TO-MAIN; some 6/6 are MAIN/SIDE-TO-SIDE
  # hence must be compatible with side chain...
  # BUT side chain is *almost* always Pro or Leu
  
  rin.B3.3 <- rin %>% 
    filter( (source_gnccn == "B3.3") ) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts, Chain.Types, file) %>% 
    unique()
  rin.B3.3 <- rin.B3.3 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
  
  
  # (3) cxb1.1-1x22
  rin.cxb1.1.1x22 <- rin %>% 
    filter( (source_gnccn == "cxb1.1") & (target_gnccn == "1x22" )) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts, Chain.Types, file) %>% 
    unique()
  rin.cxb1.1.1x22 <- rin.cxb1.1.1x22 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
  
    # 5 complexes make contacts; 5/5 are SIDE-TO-SIDE
    # hence side chains must be compatible
  
  # (4) cxb1.1-7x24
  rin.cxb1.1.7x24 <- rin %>% 
    filter( (source_gnccn == "cxb1.1") & (target_gnccn == "7x24" )) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts, Chain.Types, file) %>% 
    unique()
  rin.cxb1.1.7x24 <- rin.cxb1.1.7x24 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
  
  # 5 complexes make contacts; 5/5 are SIDE-TO-MAIN/SIDE
  # hence side chains must be compatible
 
  # (5) cxb1.1-7x24
  rin.CX.1.NTr.Cm1 <- rin %>% 
    filter( (source_gnccn == "CX.1") & (target_gnccn == "NTr.Cm1" )) %>%
    select(source_gnccn, target_gnccn, Number.of.atomic.contacts, Chain.Types, file) %>% 
    unique()
  rin.CX.1.NTr.Cm1 <- rin.CX.1.NTr.Cm1 %>% pivot_wider(names_from = Chain.Types, values_from = Number.of.atomic.contacts)
  
  