# Name:     02_motif_table.R
# Updated:  20201218
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

# **THIS SCRIPT WAS RUN IN PRIOR TO THE ORIGINAL SUBMISSION (c. 11/2020; 
# see /kleist_2020 foler) AND NOT RE-RUN DURING REVISIONS DUE TO LONG RUN-TIMES
# AND NO CHANGES TO CODE; INPUT/OUTPUT PATHS AND REQUIRED INPUTS WERE MODIFIED
# TO ACCOMODATE NEW FOLDER HIERARCHIES - ANY "DOWNSTREAM" ANALYSIS WAS RE-DONE
# AND STARTS WITH THE OUTPUT FROM THIS SCRIPT AS THE STARTING POINT

##### 1:FREQUENCY-OF-OCCURANCE TABLE ###########################################
  
  # (1) FILTERING & CLEANING
  data <- read_csv("07_ck_motif/data/processed/CK_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  
  # REMOVE DOUBLES
  data <- data %>% select(-motif_no) %>% distinct() 
  
  # ADD NUMBER OF ORTHOLOGS
  no.ortho <- read.csv("07_ck_motif/data/processed/no_ortho.csv")
  data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$ck)]
  data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$ck)]
  rm(no.ortho)
  
  # REPLACE "NA" WITH ASN-ALA
  data$motif[is.na(data$motif)] <- c("AsnAla")
  
  # REMOVE XCL1, XCL2, CCL4L1
  data <- data %>% filter(protein != "xcl1" & protein != "xcl2" & protein != "ccl4l1")
  
  # (2) CONSERVATION SCORING AT DIFFERENT LEVELS (ORTHOLOGS, FAMILY, ETC)
  # INDIVIDUAL - motif occurance across a single chemokine
  indiv <- data %>% count(motif, protein, class, mer, mask, no_ortho) %>% 
    mutate(pct_ortho = n / no_ortho)
  colnames(indiv)[6] <- c("total_ortho")
  colnames(indiv)[7] <- c("count_ortho")
  indiv <- indiv %>% select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)
  
  # FAMILY - motif occurance across each class
    # only considering motifs as occuring in a chemokine if motif appears in
    # >50% of orthologous sequences for that chemokine
    # Evaluation occurs on paralog level, not abundance of ortholog sequences
    # due to assymetric data sets for each chemokine (ie some have 50 sequences)
  family <- indiv %>% filter(pct_ortho >= 0.5) %>% count(motif, class)
  colnames(family)[3] <- c("count_family")
    
  freq <- left_join(indiv, family)
  freq[is.na(freq)] <- 0
  
  # SUPERFAMILY - motif occurance across all chemokines
  super <- indiv %>% filter(pct_ortho >= 0.5) %>% count(motif)
  colnames(super)[2] <- c("count_super")
  
  freq <- left_join(freq, super)
  freq[is.na(freq)] <- 0
  
  # ADD PCT COLUMNS
  freq <- freq %>%
    mutate(total_family = case_when(
      class == "cc" ~ 27, # recall that CCL4L1 has been removed, so 27 instead of 28
      class == "cxc" ~ 17,
      class == "cx3c" ~ 1 # note that XC family has been removed
    ))
  
  freq <- freq %>% mutate(pct_family = count_family / total_family)
  freq$total_super <- 43
  freq <- freq %>% mutate(pct_super = count_super / total_super)
  
  # (3) CLEAN AND WRITE
  freq <- freq %>% select(motif, protein, class, mer, mask, 
                          count_ortho, total_ortho, pct_ortho,
                          count_family, total_family, pct_family,
                          count_super, total_super, pct_super)
  
  write_csv(freq, "07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")
  rm(data, family, freq, indiv, super)
  