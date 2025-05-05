# Name:     02_motif_table.R
# Updated:  20210510
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

# **THIS SCRIPT WAS MODIFIED FROM THE PRIOR ANALYSIS TO COUNT MOTIF REPRESENTATION
# IN THE SUBFAMILY (CC/CXC) AND FAMILY (ALL CHEMOKINES/GPCRs) ONLY AMONG HUMAN
# PARALOGS. PREVIOUSLY COUNTS HAD BEEN DONE ACCORDING TO WHETHER A MOTIF WAS
# REPRESENTED AMONG ≥50% ORTHOLOGS REGARDLESS OF WHETHER IT WAS FOUND IN THE
# HUMAN (PARALOG) SEQUENCE. NOTE THAT THIS SCRIPT TAKES AS INPUT A FILE THAT
# WAS GENERATED DURING THE ORIGINAL SUBMISSION (c. 11/2020; see
# /kleist_2020 folder) AND NOT RE-RUN DURING REVISIONS DUE TO LONG RUN-TIMES
# AND NO CHANGES TO CODE

##### 1:FREQUENCY-OF-OCCURANCE TABLE - NTERM ###################################
  
  # setup data
  data <- read_csv("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  data <- data %>% select(-motif_no) %>% distinct() 
    # remove doubles of motifs within same sequence
  no.ortho <- read.csv("08_ckr_motif/data/processed/no_ortho.csv")
  no.ortho$protein <- as.character(no.ortho$protein)
  data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$protein)]
  data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$protein)]
  rm(no.ortho)
  data$motif[is.na(data$motif)] <- c("AsnAla") # no gaps are included so all NAs must be AsnAla
  
  # INDIVIDUAL - motif occurance across a single chemokine
  indiv <- data %>% count(motif, protein, class, mer, mask, no_ortho) %>% 
    mutate(pct_ortho = n / no_ortho)
  colnames(indiv)[6] <- c("total_ortho")
  colnames(indiv)[7] <- c("count_ortho")
  indiv <- indiv %>% select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)
  
  # Note that NA --> AsnAla and doubles already addressed
  
  # DEFINE HUMAN SEQUENCES (counts will be done among human paralogs)
  human <- data %>% filter(grepl("human", data$file))
  
  # SUBFAMILY - motif occurance across each class
    # This is changed since original submission; motifs are counted regardless
    # of ortholog conservation. Instead they are counted according to the
    # number of times they are represented in receptor subfamilies (ie CC, CXC)
    # or the receptor family (ie all chemokine receptors).
    # Evaluation occurs on paralog level, not abundance of ortholog sequences
    # due to assymetric data sets for each receptor (ie some have 50 sequences)
  family <- human %>% count(motif, class)
  colnames(family)[3] <- c("count_family")
    
  freq <- left_join(indiv, family)
  freq[is.na(freq)] <- 0
  
  
  # FAMILY - motif occurance across all receptors
  super <- human %>% count(motif)
  colnames(super)[2] <- c("count_super")
  
  freq <- left_join(freq, super)
  freq[is.na(freq)] <- 0
  
  # ADD PCT COLUMNS
  freq <- freq %>%
    mutate(total_family = case_when(
      class == "cc" ~ 10,
      class == "cxc" ~ 6,
      class == "ack" ~ 5,
      class == "xc" ~ 1,
      class == "cx3c" ~ 1
    ))
  freq <- freq %>% mutate(pct_family = count_family / total_family)
  freq$total_super <- 23
  freq <- freq %>% mutate(pct_super = count_super / total_super)
  
  freq <- freq %>% select(motif, protein, class, mer, mask, 
                          count_ortho, total_ortho, pct_ortho,
                          count_family, total_family, pct_family,
                          count_super, total_super, pct_super)
  
  # add column for whether found in human or not, filter for human only motifs
  human <- human %>% select(protein, motif)
  human$human <- c(1)
  freq <- left_join(freq, human)
  freq <- freq %>% filter(human == 1)
  freq <- freq %>% select(-human)
  
  write_csv(freq, "08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")
  rm(data, family, freq, indiv, super, human)
  
  
##### 2:FREQUENCY-OF-OCCURANCE TABLE - ECL2 ####################################
  
  # setup data
  data1 <- read_csv("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_ALL_MERS_CYSLESS.csv")
  data2 <- read_csv("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_ALL_MERS_CYSLESS.csv")
  data <- bind_rows(data1, data2)
  rm(data1, data2)
  
  # DEFINE HUMAN SEQUENCES (counts will be done among human paralogs)
  human <- data %>% filter(grepl("human", data$file))
  
  # remove doubles of motifs within same sequence
  data <- data %>% select(-file, -motif_no) %>% distinct() 
  human <- human %>% select(-file, -motif_no) %>% distinct() 
  
  # INDIVIDUAL - motif occurance across a single receptor
  # first add no ortho seqs for each receptor
  no.ortho <- read.csv("08_ckr_motif/data/processed/no_ortho.csv")
  no.ortho$protein <- as.character(no.ortho$protein)
  data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$protein)]
  data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$protein)]
  rm(no.ortho)
  data$motif[is.na(data$motif)] <- c("AsnAla") # no gaps are included so all NAs must be AsnAla
  
  # now count
  indiv <- data %>% count(motif, protein, class, mer, mask, no_ortho) %>% 
    mutate(pct_ortho = n / no_ortho)
  colnames(indiv)[6] <- c("total_ortho")
  colnames(indiv)[7] <- c("count_ortho")
  indiv <- indiv %>% select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)
  
  # FAMILY - motif occurance across each class (see above for context)
  family <- human %>% count(motif, class)
  colnames(family)[3] <- c("count_family")
  
  freq <- left_join(indiv, family)
  freq[is.na(freq)] <- 0
  
  # SUPERFAMILY - motif occurance across all chemokines
  super <- human %>% count(motif)
  colnames(super)[2] <- c("count_super")
  
  freq <- left_join(freq, super)
  freq[is.na(freq)] <- 0
  
  # ADD PCT COLUMNS
  freq <- freq %>%
    mutate(total_family = case_when(
      class == "cc" ~ 10,
      class == "cxc" ~ 6,
      class == "ack" ~ 5,
      class == "xc" ~ 1,
      class == "cx3c" ~ 1
    ))
  freq <- freq %>% mutate(pct_family = count_family / total_family)
  freq$total_super <- 23
  freq <- freq %>% mutate(pct_super = count_super / total_super)
  
  freq <- freq %>% select(motif, protein, class, mer, mask, 
                          count_ortho, total_ortho, pct_ortho,
                          count_family, total_family, pct_family,
                          count_super, total_super, pct_super)
  
  # add column for whether found in human or not, filter for human only motifs
  human <- human %>% select(protein, motif)
  human$human <- c(1)
  freq <- left_join(freq, human)
  freq <- freq %>% filter(human == 1)
  freq <- freq %>% select(-human)
  
  write_csv(freq, "08_ckr_motif/output/CKR_MOTIF_FREQUENCY_ECL2.csv")
  rm(data, family, freq, indiv, super, human)
  