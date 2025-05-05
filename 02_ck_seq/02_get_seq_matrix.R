# Name:     02_get_seq_matrix.R
# Updated:  20201020
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

# NOTE 1 (20202025):
# TOPIC: LOWER AND UPPER CASE RESIDUES
# KMAD algorithm introduced lower case letters where substitutions were made,
# however these are recognized as different symbols during classification. As such
# all symbols were made upper case before training/testing

##### FUNCTIONS ################################################################

  Align2DataFrame <- function(ALIGNMENT, NAMES){
    aln <- readAAMultipleAlignment(ALIGNMENT)
    aln.df <- as.matrix(aln)
    aln.df <- toupper(aln.df) # new 20202025
    seqname <- as.tibble(rownames(aln.df))
    colnames(seqname) <- c("seq")
    aln.df <- as.tibble(aln.df)
    gnccn <- t(read.table(NAMES, sep = ","))
    colnames(aln.df) <- gnccn
    aln.df <- cbind(seqname, aln.df)
    
    return(aln.df)
    rm(aln, aln.df, gnccn, seqname)
  }
  
##### 1: MAKE CK SEQ MATRIX LABELED WITH "CC" AND "CXC" ########################
  
  # import
  data <- Align2DataFrame("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.fasta", 
                         "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = case_when(
    grepl("cx3c" , data$protein) ~ "NA",
    grepl("cc" , data$protein) ~ "cc",
    grepl("cxc" , data$protein) ~ "cxc"
  ))
  
  # filter & reorder
  data <- data %>% filter(class == "cc" | class == "cxc")
  data <- data %>% select(protein, seq, class, NTc.Cm70:CT.308)
  
  # write output
  write_csv(data, "02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv")
  rm(data)

    
##### 2: MAKE "TEST SET" CONTAINING CHEMOKINE PARALOGS #########################
  
  # import
  data <- Align2DataFrame("02_ck_seq/data/processed/ALL_para.fasta", 
                          "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = case_when(
    grepl("cx3c" , data$protein) ~ "non",
    grepl("cc" , data$protein) ~ "cc",
    grepl("cxc" , data$protein) ~ "cxc",
    protein == "xcl1" ~ "non",
    protein == "xcl2" ~ "non"
    
  ))
  
  # reorder
  data <- data %>% select(protein, seq, class, NTc.Cm70:CT.308)
  
  # write output
  write_csv(data, "02_ck_seq/data/processed/ALL_para_df.csv")
  rm(data)
  
  
##### 3: MAKE "TEST SET" CONTAINING VIRAL CHEMOKINES ###########################
  
  # import
  data <- Align2DataFrame("02_ck_seq/data/raw/SEQUENCES_VIRAL_CK.fasta", 
                          "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, "|", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = "non")
  
  # reorder
  data <- data %>% select(protein, seq, class, NTc.Cm70:CT.308)
  
  # write
  write_csv(data, "02_ck_seq/data/processed/ALL_virus_df.csv")
  