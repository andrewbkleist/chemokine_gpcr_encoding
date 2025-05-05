# Name:     02_t2_to_bfactor.R
# Updated:  20230106
# Author:   Andrew Kleist

# packages, working directsory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)

##### 1: MAP TIER 2 SCORE TO B-FACTOR - VMIP-II ################################

  # SEQUENCE SCORES - VIRUS CC vs CXC SEQUENCE PREDICTION SCORES
  # Get tier 2 scores from viruses, select vMIP-II scores
  score <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(protein == "vmip2xhhv8p")
  # subset vMIP-II T2 scores to give only position and mean score 
  score <- score %>% select(position, mean) %>% unique()
  
  # STRUCTURE - vMIP-II:CXCR4 RESIDUE INTERFACE
  # Subset vMIP-II positions located at vMIP-II:CXCR4 interface
  # Subset vMIP-II positions that are ALSO Tier 2 (subfamily predictive)
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file == "4rws") %>%
    filter(cc_cxc_lr_ck >= 0.75) %>%
    select(source_gnccn) %>% unique()
    # get all chemokine interface positions
  rin <- rin$source_gnccn # make character vector
  
  # SEQUENCE SCORES
  # subset by positions on vMIP-II:CXCR4 interface only
  score <- score %>% filter(position %in% rin) 
  
  # (3) MAP B-FACTORS (RESIDUE SCORES) TO vMIP-II STRUCTURE
  # import PDB, add column to PDB reflecting CCN positions
  pdb.file <- read.pdb("35_viral_ackr/data/processed/vmipii.pdb")
  pdb <- pdb.file$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>%
    select(ccn_4rws_ck, clean_4rws_ck)
  colnames(lookup) <- c("position", "resno")

  # map CCN positions to pdb residue numbers
  pdb <- left_join(pdb, lookup)
  rm(lookup)
  
  # Map Prediction Probability Scores (PPS) to B-factor column in PDB
  # Using CCN positions as reference/lookup table
  # Note: the actual numbers that are mapped as B-factors (PPS) 
  # not the "Tier 2" scores themselves (which is a predictive accuracy score)
  # but the score given from the logistic regression model, wherein scores
  # close to 1 indicate a prediction of CXC-like residue and scores close to
  # 0 indicate prediction of a CC-like residue
  pdb <- left_join(pdb, score)
  pdb$b <- pdb$mean
  pdb$b <- round(pdb$b, 2)
  
  # summary
  # test <- pdb %>% select(b, position) %>% unique()

  # Provide B-factors for non-interface, non-Tier 2 residues
  # Arbitrarily setting PPS to 0.5 since this would indicate 
  # no preference CC vs CXC
  pdb <- pdb %>% mutate(b = case_when(
    is.na(b) ~ 0.5,
    !is.na(b) ~b
  ))
  
  # remove CCN position from PDB file
  pdb <- pdb %>% select(-position, - mean)
  
  # write pdb
  pdb.file$atom <- pdb
  # write.pdb(pdb.file, file = "35_viral_ackr/output/vmipii_t2_score_bfactor_20230106.pdb")
  
  
##### 2: MAP TIER 2 SCORE TO B-FACTOR - CCL5 ###################################
  
  # SEQUENCE SCORES - VIRUS CC vs CXC SEQUENCE PREDICTION SCORES
  # Get tier 2 scores, select CCL5 scores
  score <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(protein == "ccl5")
  # subset CCL5 T2 scores to give only position and mean score 
  score <- score %>% select(position, mean) %>% unique()
  
  # STRUCTURE - CCL5:CCR5 RESIDUE INTERFACE
  # Subset vMIP-II positions located at CCL5:CCR5 (zheng) interface
  # Subset CCL5 positions that are ALSO Tier 2 (subfamily predictive)
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file == "zheng") %>%
    filter(cc_cxc_lr_ck >= 0.75) %>%
    select(source_gnccn) %>% unique()
  # get all chemokine interface positions
  rin <- rin$source_gnccn # make character vector
  
  # SEQUENCE SCORES
  # subset by positions on vMIP-II:CXCR4 interface only
  score <- score %>% filter(position %in% rin) 
  
  # (3) MAP B-FACTORS (RESIDUE SCORES) TO CCL5 STRUCTURE
  # import PDB, add column to PDB reflecting CCN positions
  pdb.file <- read.pdb("35_viral_ackr/data/processed/ccl5.pdb")
  pdb <- pdb.file$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>%
    select(ccn_zheng_ck, clean_zheng_ck)
  colnames(lookup) <- c("position", "resno")
  
  # map CCN positions to pdb residue numbers
  pdb <- left_join(pdb, lookup)
  rm(lookup)
  
  # Map Prediction Probability Scores (PPS) to B-factor column in PDB
  # Using CCN positions as reference/lookup table
  # Note: the actual numbers that are mapped as B-factors (PPS) 
  # not the "Tier 2" scores themselves (which is a predictive accuracy score)
  # but the score given from the logistic regression model, wherein scores
  # close to 1 indicate a prediction of CXC-like residue and scores close to
  # 0 indicate prediction of a CC-like residue
  pdb <- left_join(pdb, score)
  pdb$b <- pdb$mean
  pdb$b <- round(pdb$b, 2)
  
  # summary
  # test <- pdb %>% select(b, position) %>% unique()
  
  # Provide B-factors for non-interface, non-Tier 2 residues
  # Arbitrarily setting PPS to 0.5 since this would indicate 
  # no preference CC vs CXC
  pdb <- pdb %>% mutate(b = case_when(
    is.na(b) ~ 0.5,
    !is.na(b) ~b
  ))
  
  # remove CCN position from PDB file
  pdb <- pdb %>% select(-position, - mean)
  
  # write pdb
  pdb.file$atom <- pdb
  # write.pdb(pdb.file, file = "35_viral_ackr/output/ccl5_t2_score_bfactor_20230106.pdb")
  