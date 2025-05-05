# Name:     01_map_satmut_to_rins.R
# Updated:  20210514
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)
library(bio3d)

##### 1: MAP SATMUT TO PDB #####################################################

  # import rins, import satmut and average, combine with rins
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") 
  rin <- rin %>% filter(file == "ngo") 
  satmut <- read_csv("09_sat_mut_data/output/cxcr4_clean_gpcrdb_means_log.csv") %>% 
    select(gn, sele, sub_mean) %>%
    filter(sele == "cxcl12") %>%
    group_by(gn) %>%
    summarise(sub_mean2 = mean(sub_mean)) %>%
    ungroup()
  colnames(satmut)[1] <- c("target_gnccn")
  rin <- left_join(rin, satmut)

  # map gnccn to pdb
  pdb <- read.pdb("01_structure_contacts/data/pdbs/ngo_model_cxcr4_clean.pdb", multi=TRUE)
  pdb.atom <- pdb$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>% select(bw_ngo_ckr, clean_ngo_ckr)
  pdb.atom$target_gnccn <- lookup$bw_ngo_ckr[match(unlist(pdb.atom$resno), lookup$clean_ngo_ckr)]
  
  # map values to pdb b factor
  satmut <- satmut %>% select(target_gnccn, sub_mean2) %>% unique()
  pdb.atom$b <- satmut$sub_mean2[match(unlist(pdb.atom$target_gnccn), satmut$target_gnccn)]
  pdb.atom$b <- abs(pdb.atom$b)
  pdb.atom <- pdb.atom %>% mutate(b = case_when(
    is.na(b) ~ 0,
    !is.na(b) ~ b
  ))
  
  # write pdb
  write.pdb(pdb=pdb, b = pdb.atom$b, file = "54_cxcr4_sat_mut/output/cxcr4_satmut_bf.pdb")
  
  