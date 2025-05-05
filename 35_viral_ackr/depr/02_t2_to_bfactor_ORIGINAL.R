# Name:     02_t2_to_bfactor.R
# Updated:  20230106
# Author:   Andrew Kleist

# packages, working directsory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)

##### 1: MAP TIER 2 SCORE TO B-FACTOR - VMIP-II ################################

  # SCOREFILE, PROTEIN, RINFILE, 
  # Get T2 scores from viruses, select vMIP-II
  score <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(protein == "vmip2xhhv8p")
  
  # Get T2 scores all positions
  t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean >= 0.75)
    # T2 scores defined as mean score >= 0.75
  
  # Get T2 at interface
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    select(source_gnccn) %>% unique()
    # get all chemokine interface positions
  rin <- rin$source_gnccn
  t2 <- t2 %>% filter(motif %in% rin) 
    # filter T2 to select only those occuring at interface
  rm(rin)
  t2 <- t2$motif
  
  # Get T2 scores located at interface residues
  score <- score %>% filter(position %in% t2)
  rm(t2)
  
  # subset vMIP-II T2 scores to give only position and mean score 
  score <- score %>% select(position, mean) %>% unique()
  
  # import, convert PDB
  pdb.file <- read.pdb("35_viral_ackr/data/processed/vmipii.pdb")
  pdb <- pdb.file$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>%
    select(ccn_4rws_ck, clean_4rws_ck)
  colnames(lookup) <- c("position", "resno")

  # map CCN positions to pdb residue numbers
  pdb <- left_join(pdb, lookup)
  rm(lookup)
  
  # Map bfactor to PDB residue numbers using CCN positions as lookup table;
  # Note that  B-factors will only be mapped for Tier 2 residues,
  # ie residues for which T2 scores are >= 0.75.
  # Importantly, the actual numbers that are mapped as B-factors are 
  # not the T2 scores themselves (which is a predictive accuracy score)
  # but the prediction probbability score from the logistic regression model.
  # These scores give a prediction, for every 
  pdb <- left_join(pdb, score)
  pdb$b <- pdb$mean
  pdb$b <- round(pdb$b, 2)

  pdb <- pdb %>% select(-position, - mean)    
  pdb <- pdb %>% mutate(b = case_when(
    is.na(b) ~ 0.5,
    !is.na(b) ~b
  ))
  
  # write pdb
  pdb.file$atom <- pdb
  # write.pdb(pdb.file, file = "35_viral_ackr/output/vmipii_t2_score_bfactor.pdb")
  
  
##### 2: MAP TIER 2 SCORE TO B-FACTOR - CCL5 ###################################
  
  # Get T2 scores, select protein
  score <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(protein == "ccl5")
  
  # Get T2
  t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean > 0.75)
  
  # Get T2 at interface
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    select(source_gnccn) %>% unique()
  rin <- rin$source_gnccn
  t2 <- t2 %>% filter(motif %in% rin)
  rm(rin)
  t2 <- t2$motif
  
  # Get score at T2 and interface
  score <- score %>% filter(position %in% t2)
  rm(t2)
  
  # simplify 
  score <- score %>% select(position, mean) %>% unique()
  
  # import, convert PDB
  pdb.file <- read.pdb("35_viral_ackr/data/processed/ccl5.pdb")
  pdb <- pdb.file$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>%
    select(ccn_zheng_ck, clean_zheng_ck)
  colnames(lookup) <- c("position", "resno")
  
  # map to pdb
  pdb <- left_join(pdb, lookup)
  rm(lookup)
  
  # map bfactor
  pdb <- left_join(pdb, score)
  pdb$b <- pdb$mean
  pdb$b <- round(pdb$b, 2)
  
  pdb <- pdb %>% select(-position, - mean)    
  pdb <- pdb %>% mutate(b = case_when(
    is.na(b) ~ 0.5,
    !is.na(b) ~b
  ))
  
  # write pdb
  pdb.file$atom <- pdb
  write.pdb(pdb.file, file = "35_viral_ackr/output/ccl5_t2_score_bfactor.pdb")
  
  
##### 1: MAP TIER 2 SCORE TO B-FACTOR - VMIP-II ################################
  
  # SCOREFILE, PROTEIN, RINFILE, 
  # Get T2 scores from viruses, select vMIP-II
  score <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(protein == "vmip2xhhv8p")
  
  # subset vMIP-II T2 scores to give only position and mean score 
  score <- score %>% select(position, mean) %>% unique()
  
  # import, convert PDB
  pdb.file <- read.pdb("35_viral_ackr/data/processed/vmipii.pdb")
  pdb <- pdb.file$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20220729.csv") %>%
    select(ccn_4rws_ck, clean_4rws_ck)
  colnames(lookup) <- c("position", "resno")
  
  # map CCN positions to pdb residue numbers
  pdb <- left_join(pdb, lookup)
  rm(lookup)
  
  # Map bfactor to PDB residue numbers using CCN positions as lookup table
  pdb <- left_join(pdb, score)
  pdb$b <- pdb$mean
  pdb$b <- round(pdb$b, 2)
  
  pdb <- pdb %>% select(-position, - mean)    
  # pdb <- pdb %>% mutate(b = case_when(
  #   is.na(b) ~ 0.5,
  #   !is.na(b) ~b
  # ))
  # 
  # write pdb
  pdb.file$atom <- pdb
  write.pdb(pdb.file, file = "35_viral_ackr/output/vmipii_t2_score_bfactor_TEST.pdb")
  
  
  
  
  
  
  
  
  
  
  
  ##### 2: MAP TIER 2 SCORE TO B-FACTOR - CCL5 ###################################
  
  # Get T2 scores, select protein
  score <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(protein == "ccl5")
  
  # Get T2
  t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
    filter(mean > 0.75)
  
  # Get T2 at interface
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    select(source_gnccn) %>% unique()
  rin <- rin$source_gnccn
  t2 <- t2 %>% filter(motif %in% rin)
  rm(rin)
  t2 <- t2$motif
  
  # Get score at T2 and interface
  score <- score %>% filter(position %in% t2)
  rm(t2)
  
  # simplify 
  score <- score %>% select(position, mean) %>% unique()
  
  # import, convert PDB
  pdb.file <- read.pdb("35_viral_ackr/data/processed/ccl5.pdb")
  pdb <- pdb.file$atom
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>%
    select(ccn_zheng_ck, clean_zheng_ck)
  colnames(lookup) <- c("position", "resno")
  
  # map to pdb
  pdb <- left_join(pdb, lookup)
  rm(lookup)
  
  # map bfactor
  pdb <- left_join(pdb, score)
  pdb$b <- pdb$mean
  pdb$b <- round(pdb$b, 2)
  
  pdb <- pdb %>% select(-position, - mean)    
  pdb <- pdb %>% mutate(b = case_when(
    is.na(b) ~ 0.5,
    !is.na(b) ~b
  ))
  
  # write pdb
  pdb.file$atom <- pdb
  write.pdb(pdb.file, file = "35_viral_ackr/output/ccl5_t2_score_bfactor.pdb")
  
  
  
  
  
  
  
  
  
  