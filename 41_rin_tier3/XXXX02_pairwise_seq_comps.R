# Name:     02_pairwise_seq_comps.R
# Updated:  20210408
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)

##### FUNCTIONS ################################################################

  GetPairwiseEntropy <- function(SEQS, ALN, NAME){
    # import aln, convert to matrix
    seqs <- c(SEQS)
    aln <- read.fasta(ALN)
    aln <- as.data.frame(aln$ali)
    aln$seq <- rownames(aln)
    aln <- aln %>% filter(seq %in% seqs)
    aln <- aln %>% select(-seq)
    
    # convert to matrix and score
    aln <- as.matrix(aln)
    # test <- conserv(x=aln, method="similarity", sub.matrix="bio3d")  
    aln   <- entropy(aln)
    aln <- as.data.frame(aln$H.10.norm)
    colnames(aln) <- c(NAME)
    
    # add common numbering
    # assign column names
    aln.names <- c(read.table("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
    aln.names <- as.data.frame(aln.names)
    aln.names <- t(aln.names)
    aln.names <- as.data.frame(aln.names)
    colnames(aln.names) <- c("ccn")
    aln <- bind_cols(aln.names, aln)
    rownames(aln) <- 1:nrow(aln)
    rm(aln.names)
    return(aln)
    rm(aln)
  }
  
  
##### 1: GET PAIRWISE COMPARISONS ##############################################
  
  # 1
  ccl5.ccl28 <- GetPairwiseEntropy(c("ccl5_human", "ccl28_human"), 
                                   "02_ck_seq/data/processed/ALL_para_simple_name.fasta",
                                   "entropy10")
  # ccl5.ccl28 <- ccl5.ccl28 %>% unite(col = "ccn_01", ccn:entropy10, sep = "_")
  # ccl5.ccl28$entropy10 <- as.character(ccl5.ccl28$entropy10)
  
  # 2
  ccl5.ccl11 <- GetPairwiseEntropy(c("ccl5_human", "ccl11_human"), 
                                   "02_ck_seq/data/processed/ALL_para_simple_name.fasta",
                                   "entropy10")
  # ccl5.ccl11 <- ccl5.ccl11 %>% unite(col = "ccn_01", ccn:entropy10, sep = "_")
  
  # ccl5.ccl11$entropy10 <- as.character(ccl5.ccl11$entropy10)
  
  
  # 3
  ccl27.ccl28 <- GetPairwiseEntropy(c("ccl27_human", "ccl28_human"), 
                                   "02_ck_seq/data/processed/ALL_para_simple_name.fasta",
                                   "entropy10")
  # ccl27.ccl28$entropy10 <- as.character(ccl27.ccl28$entropy10)
  # ccl27.ccl28 <- ccl27.ccl28 %>% unite(col = "ccn_01", ccn:entropy10, sep = "_")
  
  
  
  # combine
  data <- data.frame(ccn = ccl5.ccl11$ccn, 
                     ccl5_ccl11 = ccl5.ccl11$entropy10, 
                     ccl5_ccl28 = ccl5.ccl28$entropy10,
                     ccl27_ccl28 = ccl27.ccl28$entropy10)
  
  # 
  ccr10 <- data %>% filter((ccl27_ccl28 == 1) & (ccl5_ccl28 == 0) )
  ccr3 <- data %>% filter((ccl27_ccl28 == 0) & (ccl5_ccl28 == 1) )
  
  
  # setdiff etc
  test <- setdiff(ccl5.ccl28, ccl27.ccl28)
  
  
  
  
  
  

  
  # test <- as.data.frame(aln$ali[,100]) # XCL1, 46th resi column
  # test <- as.data.frame(aa123( aln$ali[c( "ccl5_human/1-91", "xcl2_human/1-114"),] ))
  
  
  
  test   <- entropy(aln)
  test <- as.data.frame(test$H.10.norm)
  
  test <- conserv(x=aln, method="similarity", sub.matrix="bio3d")  
  
  
  
  
  test   <- entropy(aln)
  