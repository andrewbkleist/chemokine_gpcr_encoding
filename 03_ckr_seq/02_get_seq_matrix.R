# Name:     02_get_seq_matrix.R
# Updated:  20210330
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### FUNCTIONS ################################################################

  Align2DataFrame <- function(ALIGNMENT, NAMES){
    aln <- readAAMultipleAlignment(ALIGNMENT)
    aln.df <- as.matrix(aln)
    seqname <- as.tibble(rownames(aln.df))
    colnames(seqname) <- c("seq")
    aln.df <- as.tibble(aln.df)
    gnccn <- t(read.table(NAMES, sep = ","))
    colnames(aln.df) <- gnccn
    aln.df <- cbind(seqname, aln.df)
    
    return(aln.df)
    rm(aln, aln.df, gnccn, seqname)
  }

##### 1: MAKE CKR SEQ MATRIX LABELED WITH "CC" AND "CXC" #######################
  
  # import
  data <- Align2DataFrame("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta", 
                         "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = case_when(
    grepl("ccrl2" , data$protein) ~ "NA",
    grepl("cc" , data$protein) ~ "cc",
    grepl("cxc" , data$protein) ~ "cxc"
  ))
  
  # filter & reorder
  data <- data %>% filter(class == "cc" | class == "cxc")
  data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)
  
  # write output
  write_csv(data, "03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv")
  rm(data)
  
  
##### 2: MAKE NEW "TEST SET" CONTAINING CLASS A GPCRS ##########################
  
  # import
  data <- Align2DataFrame("03_ckr_seq/data/processed/ALL_classa.fasta", 
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, "_", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = case_when(
    grepl("ccrl2" , data$protein) ~ "non",
    grepl("cckar" , data$protein) ~ "non",
    grepl("cc" , data$protein) ~ "cc",
    grepl("cxc" , data$protein) ~ "cxc"
  ))
  
  data <- data %>% mutate(class = case_when(
    is.na(class) ~ "non",
    !is.na(class) ~ class
  ))
  
  # reorder
  data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)
  
  # write
  write_csv(data, "03_ckr_seq/data/processed/ALL_classa_df.csv")
  
  
##### 3: MAKE 2nd "TEST SET" CONTAINING VIRAL GPCRS ############################
  
  # import
  data <- Align2DataFrame("03_ckr_seq/data/processed/ALL_virus.fasta", 
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, "|", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = "non")
  
  # reorder
  data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)
  
  # write
  write_csv(data, "03_ckr_seq/data/processed/ALL_virus_df.csv")
  
  
##### 4: MAKE SET CONTAIING ALL PARALOGS MINUS ACKRs ###########################
  
  # import
  data <- Align2DataFrame("03_ckr_seq/data/processed/ALL_non_ackr_para.fasta", 
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")
  
  # add protein name
  pname <- data %>% separate(seq, "|", remove = TRUE) %>% select(1) 
  colnames(pname) <- c("protein")
  data <- cbind(pname, data)
  rm(pname)
  
  # cc/cxc classifier
  data <- data %>% mutate(class = "non")
  
  # reorder
  data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)
  
  # write
  write_csv(data, "03_ckr_seq/data/processed/ALL_non_ackr_para_df.csv")

  
# ##### 4: MAKE MATRIX WITH ***CCL7 BINDING RECEPTORS***  ########################
# # 285 non-CCL7-binding sequences /  136 CCL7-binding sequences 
#   # import
#   data <- Align2DataFrame("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta", 
#                           "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")
#   
#   # add protein name
#   pname <- data %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
#   colnames(pname) <- c("protein")
#   data <- cbind(pname, data)
#   rm(pname)
#   
#   # select/remove
#   data <- data %>% filter(grepl("ccr" , data$protein))
#   data <- data %>% filter(!grepl("ccrl2" , data$protein))
#   
#   
#   # classifier
#   data <- data %>% mutate(class = case_when(
#     protein == "ccr1" ~ "ccl7",
#     protein == "ccr2" ~ "ccl7",
#     protein == "ccr3" ~ "ccl7",
#     protein == "ccr5" ~ "ccl7"
#   ))
#   
#   data <- data %>% mutate(class = case_when(
#     !is.na(class) ~ class,
#     is.na(class) ~ "non_ccl7"
#   ))
#   
#   # filter & reorder
#   data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)
#   
#   # write output
#   write_csv(data, "03_ckr_seq/data/processed/ALL_ccl7_df.csv")
#   rm(data)
#   
#   
# ##### 5: MAKE MATRIX WITH ***CXCL11 BINDING RECEPTORS***  ######################
# # 682 non-CXCL11-binding sequences /  175 CCL7-binding sequences 
# 
#   # import
#   data <- Align2DataFrame("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta", 
#                           "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")
#   
#   # add protein name
#   pname <- data %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
#   colnames(pname) <- c("protein")
#   data <- cbind(pname, data)
#   rm(pname)
#   
#   # select/remove
#   data <- data %>% filter(!grepl("cx3c" , data$protein))
#   data <- data %>% filter(!grepl("xxc" , data$protein))
# 
#   # classifier
#   data <- data %>% mutate(class = case_when(
#     protein == "ccr3" ~ "cxcl11",
#     protein == "ccr5" ~ "cxcl11",
#     protein == "cxcr3" ~ "cxcl11",
#     protein == "ackr1" ~ "cxcl11",
#     protein == "ackr2" ~ "cxcl11"
#   ))
#   
#   data <- data %>% mutate(class = case_when(
#     !is.na(class) ~ class,
#     is.na(class) ~ "non_cxcl11"
#   ))
#   
#   # filter & reorder
#   data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)
#   
#   # write output
#   write_csv(data, "03_ckr_seq/data/processed/ALL_cxcl11_df.csv")
#   rm(data)
