# Name:     04_cons_seq_logo.R
# Updated:  20210416
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggseqlogo)
require(Biostrings)


##### FUNCTIONS  ###############################################################

  GetLogoPosition <- function(ALIGNMENT, LOOKUP, POSITION){
    aln <- readAAMultipleAlignment(ALIGNMENT)
    aln.df <- as.matrix(aln)
    aln.df <- toupper(aln.df) # new 20202025
    seqname <- as.tibble(rownames(aln.df))
    colnames(seqname) <- c("seq")
    aln.df <- as.tibble(aln.df)
    gnccn <- t(read.table(LOOKUP, sep = ","))
    colnames(aln.df) <- gnccn
    aln.df <- cbind(seqname, aln.df)

    rm(aln, gnccn, seqname)
    
    # manually define positions for seq logo
    pos <- c(POSITION)
    
    # load alignment df
    # aln.df <- read_csv("02_ck_seq/data/processed/ALL_para_df.csv")
    aln.df <- aln.df[(names(aln.df) %in% pos)]
    # aln.df <- unite(aln.df[,1:ncol(aln.df)], col = seq_string,  sep = "")
    return(aln.df)
  }
  
##### 1: GET LOGOS #############################################################
  
  # (1.1) CHEMOKINE ------------------------------------------------------------
  
  # t1
  logo <- GetLogoPosition("02_ck_seq/data/processed/ALL_para.fasta",
                          "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt",
                          "b1b2.9")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  # t2
  logo <- GetLogoPosition("02_ck_seq/data/processed/CC_para.fasta",
                          "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt",
                          "b1b2.9")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  # t3
  logo <- GetLogoPosition("02_ck_seq/data/processed/ccl5_ortho.fasta",
                          "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt",
                          "b1b2.9")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  # t3
  logo <- GetLogoPosition("02_ck_seq/data/processed/ccl20_ortho.fasta",
                          "02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt",
                          "b1b2.9")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  
  # (1.2) RECEPTOR -------------------------------------------------------------
  
  # t1
  logo <- GetLogoPosition("03_ckr_seq/data/processed/ALL_non_ackr_para.fasta",
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt",
                          "gnECL2.Cp4")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  # t1
  logo <- GetLogoPosition("03_ckr_seq/data/processed/ALL_non_ackr_para.fasta",
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt",
                          "gnECL2.Cp3")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  # t2
  logo <- GetLogoPosition("03_ckr_seq/data/processed/CC_para.fasta",
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt",
                          "gnECL2.Cp4")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  # t3
  logo <- GetLogoPosition("03_ckr_seq/data/processed/ccr5_ortho.fasta",
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt",
                          "gnECL2.Cp4")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  
  # t3
  logo <- GetLogoPosition("03_ckr_seq/data/processed/ccr6_ortho.fasta",
                          "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt",
                          "gnECL2.Cp3")
  logo <- ggseqlogo(logo,  method = 'bits' ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  gridExtra::grid.arrange(logo)
  
  
  
 
  
  