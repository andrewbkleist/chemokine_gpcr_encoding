# Name:     01_get_conservation_scores.R
# Updated:  20210323
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

# NOTE 1 (20200924):
# N-TERMINUS & CCN 
# Selected NT.1:NT.95 from master alignment of chemokine orthologs (~1056)
# and wrote to a fasta. Removed gaps relative to N-terminal Cys but did so by 
# being consistent across all orthologs for single chemokine ("column wise removal")
# Saved and re-imported, combined with core alignment to make new master align
# Saved as SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.fasta in
# chemokine_gpcr_ms/02_ck_seq/data/raw/
# Altered CCN Nterm numbering to reflect NTc.Cm70-NTc.Cm1 reflecting that
# numberung is now only relative to Cys
#
#
# NOTE 2 (20200924):
# TOPIC: ALIGNMENT SETS
# Want the following per-position descriptors:
# - Ortholog conservation (x43)
# - CC paralog conservation (1)
# - CXC paralog conservation (1)
# - CC and CXC paralog conservation (1)
# - All human (46 seq) paralog conservation (1)
# - All 1056 paralog conservation (1) - REMOVED, 20210303 NOTE
# - T2 scores (1) [different script]
# - Chemokine-specific T2 scores (viral chemokines, x ??) [different script]
# **For alignment sets - manually created the different sequence sets
#
#
# NOTE 3 (20201022):
# TOPIC: VIRAL CHEMOKINES
# (1) Sequences from Alex Gunnarsson; includes 30 sequences
# (2) Combined 30 viral sequences with 46 human paralogs and aligned using 
#     ClustalO online: viral_human_ckrs_clustal.fasta
# (3) Cross checked following alignment pairs by comparing to structural
#     alignmnet: vMIP2-CCL2 (good)
# (4) Pasted 30 viral sequences plus CCL2 sequence into new alignment with
#     46 human paralogs from master alignment; used CCL2 as bridge to adjust
#     viral sequences to master alignment; ***insertions in viral sequences
#     relative to master alignment were deleted***
# (5) Final alignment: SEQUENCES_VIRAL_CK.fasta and placed in 
#     /02_ck_seq/data/raw
#
#
# NOTE 4 (20210303):
# TOPIC: ERROR WITH MSTATX CALC > 500 SEQ
# - Noted that MstatX has 500 seq maximum
# - No other calculations should not be affected (ortholog conservation score is
#   the only one using >500 seq)
# - Addendum (20230817): use of Mstatx for >500 sequence was remedied
#
##### PART 1: CONSERVATION SCORING #############################################
# From all alignments, calculate per-position conservation scores
# Note that all sequence sets were manually curated
# the following commands were run in the command line from:
# /Users/akleist/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/02_ck_seq/data/processed
#
# mstatx -i ALL_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_para.txt
# REMOVED: mstatx -i ALL_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_ortho.txt
# mstatx -i CC_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o CC_para.txt
# mstatx -i CXC_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o CXC_para.txt
# mstatx -i ALL_CC_CXC_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_CC_CXC_para.txt
# REMOVED: mstatx -i ALL_CC_CXC_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_CC_CXC_ortho.txt
# mstatx -i ccl1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl1_ortho.txt
# mstatx -i ccl2_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl2_ortho.txt
# mstatx -i ccl3_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl3_ortho.txt
# mstatx -i ccl3l1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl3l1_ortho.txt
# mstatx -i ccl4_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl4_ortho.txt
# mstatx -i ccl5_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl5_ortho.txt
# mstatx -i ccl7_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl7_ortho.txt
# mstatx -i ccl8_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl8_ortho.txt
# mstatx -i ccl11_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl11_ortho.txt
# mstatx -i ccl13_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl13_ortho.txt
# mstatx -i ccl14_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl14_ortho.txt
# mstatx -i ccl15_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl15_ortho.txt
# mstatx -i ccl16_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl16_ortho.txt
# mstatx -i ccl17_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl17_ortho.txt
# mstatx -i ccl18_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl18_ortho.txt
# mstatx -i ccl19_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl19_ortho.txt
# mstatx -i ccl20_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl20_ortho.txt
# mstatx -i ccl21_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl21_ortho.txt
# mstatx -i ccl22_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl22_ortho.txt
# mstatx -i ccl23_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl23_ortho.txt
# mstatx -i ccl24_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl24_ortho.txt
# mstatx -i ccl25_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl25_ortho.txt
# mstatx -i ccl26_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl26_ortho.txt
# mstatx -i ccl27_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl27_ortho.txt
# mstatx -i ccl28_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccl28_ortho.txt
# mstatx -i cxcl1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl1_ortho.txt
# mstatx -i cxcl2_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl2_ortho.txt
# mstatx -i cxcl3_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl3_ortho.txt
# mstatx -i cxcl4_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl4_ortho.txt
# mstatx -i cxcl4l1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl4l1_ortho.txt
# mstatx -i cxcl5_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl5_ortho.txt
# mstatx -i cxcl6_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl6_ortho.txt
# mstatx -i cxcl7_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl7_ortho.txt
# mstatx -i cxcl8_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl8_ortho.txt
# mstatx -i cxcl9_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl9_ortho.txt
# mstatx -i cxcl10_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl10_ortho.txt
# mstatx -i cxcl11_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl11_ortho.txt
# mstatx -i cxcl12_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl12_ortho.txt
# mstatx -i cxcl13_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl13_ortho.txt
# mstatx -i cxcl14_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl14_ortho.txt
# mstatx -i cxcl16_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl16_ortho.txt
# mstatx -i cxcl17_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcl17_ortho.txt
# mstatx -i cx3cl1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cx3cl1_ortho.txt


##### PART 2: IMPORT AND TIDY SCORES ###########################################
  
  # (2.1) IMPORT AND TIDY ------------------------------------------------------
  # import scores
  all.para <- read.table("02_ck_seq/data/processed/ALL_para.txt", header = FALSE)
  all.para <- all.para %>% select(-V1)
  colnames(all.para) <- c("all_para")
  
  all.cc.cxc.para <- read.table("02_ck_seq/data/processed/ALL_CC_CXC_para.txt", header = FALSE)
  all.cc.cxc.para <- all.cc.cxc.para %>% select(-V1)
  colnames(all.cc.cxc.para) <- c("all_cc_cxc_para")
  
  cc.para <- read.table("02_ck_seq/data/processed/CC_para.txt", header = FALSE)
  cc.para <- cc.para %>% select(-V1)
  colnames(cc.para) <- c("cc_para")
  
  cxc.para <- read.table("02_ck_seq/data/processed/CXC_para.txt", header = FALSE)
  cxc.para <- cxc.para %>% select(-V1)
  colnames(cxc.para) <- c("cxc_para")
  
  # REMOVED
  # all.ortho <- read.table("02_ck_seq/data/processed/ALL_ortho.txt", header = FALSE)
  # all.ortho <- all.ortho %>% select(-V1)
  # colnames(all.ortho) <- c("all_ortho")
  
  # REMOVED
  # all.cc.cxc.ortho <- read.table("02_ck_seq/data/processed/ALL_CC_CXC_ortho.txt", header = FALSE)
  # all.cc.cxc.ortho <- all.cc.cxc.ortho %>% select(-V1)
  # colnames(all.cc.cxc.ortho) <- c("all_cc_cxc_ortho")
  
  ccl1.ortho <- read.table("02_ck_seq/data/processed/ccl1_ortho.txt", header = FALSE)
  ccl1.ortho <- ccl1.ortho %>% select(-V1)
  colnames(ccl1.ortho) <- c("ccl1_ortho")
  
  ccl2.ortho <- read.table("02_ck_seq/data/processed/ccl2_ortho.txt", header = FALSE)
  ccl2.ortho <- ccl2.ortho %>% select(-V1)
  colnames(ccl2.ortho) <- c("ccl2_ortho")
  
  ccl3.ortho <- read.table("02_ck_seq/data/processed/ccl3_ortho.txt", header = FALSE)
  ccl3.ortho <- ccl3.ortho %>% select(-V1)
  colnames(ccl3.ortho) <- c("ccl3_ortho")
  
  ccl3l1.ortho <- read.table("02_ck_seq/data/processed/ccl3l1_ortho.txt", header = FALSE)
  ccl3l1.ortho <- ccl3l1.ortho %>% select(-V1)
  colnames(ccl3l1.ortho) <- c("ccl3l1_ortho")
  
  ccl4.ortho <- read.table("02_ck_seq/data/processed/ccl4_ortho.txt", header = FALSE)
  ccl4.ortho <- ccl4.ortho %>% select(-V1)
  colnames(ccl4.ortho) <- c("ccl4_ortho")
  
  ccl5.ortho <- read.table("02_ck_seq/data/processed/ccl5_ortho.txt", header = FALSE)
  ccl5.ortho <- ccl5.ortho %>% select(-V1)
  colnames(ccl5.ortho) <- c("ccl5_ortho")
  
  ccl7.ortho <- read.table("02_ck_seq/data/processed/ccl7_ortho.txt", header = FALSE)
  ccl7.ortho <- ccl7.ortho %>% select(-V1)
  colnames(ccl7.ortho) <- c("ccl7_ortho")
  
  ccl8.ortho <- read.table("02_ck_seq/data/processed/ccl8_ortho.txt", header = FALSE)
  ccl8.ortho <- ccl8.ortho %>% select(-V1)
  colnames(ccl8.ortho) <- c("ccl8_ortho")
  
  ccl11.ortho <- read.table("02_ck_seq/data/processed/ccl11_ortho.txt", header = FALSE)
  ccl11.ortho <- ccl11.ortho %>% select(-V1)
  colnames(ccl11.ortho) <- c("ccl11_ortho")
  
  ccl13.ortho <- read.table("02_ck_seq/data/processed/ccl13_ortho.txt", header = FALSE)
  ccl13.ortho <- ccl13.ortho %>% select(-V1)
  colnames(ccl13.ortho) <- c("ccl13_ortho")
  
  ccl14.ortho <- read.table("02_ck_seq/data/processed/ccl14_ortho.txt", header = FALSE)
  ccl14.ortho <- ccl14.ortho %>% select(-V1)
  colnames(ccl14.ortho) <- c("ccl14_ortho")
  
  ccl15.ortho <- read.table("02_ck_seq/data/processed/ccl15_ortho.txt", header = FALSE)
  ccl15.ortho <- ccl15.ortho %>% select(-V1)
  colnames(ccl15.ortho) <- c("ccl15_ortho")
  
  ccl16.ortho <- read.table("02_ck_seq/data/processed/ccl16_ortho.txt", header = FALSE)
  ccl16.ortho <- ccl16.ortho %>% select(-V1)
  colnames(ccl16.ortho) <- c("ccl16_ortho")
  
  ccl17.ortho <- read.table("02_ck_seq/data/processed/ccl17_ortho.txt", header = FALSE)
  ccl17.ortho <- ccl17.ortho %>% select(-V1)
  colnames(ccl17.ortho) <- c("ccl17_ortho")
  
  ccl18.ortho <- read.table("02_ck_seq/data/processed/ccl18_ortho.txt", header = FALSE)
  ccl18.ortho <- ccl18.ortho %>% select(-V1)
  colnames(ccl18.ortho) <- c("ccl18_ortho")
  
  ccl19.ortho <- read.table("02_ck_seq/data/processed/ccl19_ortho.txt", header = FALSE)
  ccl19.ortho <- ccl19.ortho %>% select(-V1)
  colnames(ccl19.ortho) <- c("ccl19_ortho")
  
  ccl20.ortho <- read.table("02_ck_seq/data/processed/ccl20_ortho.txt", header = FALSE)
  ccl20.ortho <- ccl20.ortho %>% select(-V1)
  colnames(ccl20.ortho) <- c("ccl20_ortho")
  
  ccl21.ortho <- read.table("02_ck_seq/data/processed/ccl21_ortho.txt", header = FALSE)
  ccl21.ortho <- ccl21.ortho %>% select(-V1)
  colnames(ccl21.ortho) <- c("ccl21_ortho")
  
  ccl22.ortho <- read.table("02_ck_seq/data/processed/ccl22_ortho.txt", header = FALSE)
  ccl22.ortho <- ccl22.ortho %>% select(-V1)
  colnames(ccl22.ortho) <- c("ccl22_ortho")
  
  ccl23.ortho <- read.table("02_ck_seq/data/processed/ccl23_ortho.txt", header = FALSE)
  ccl23.ortho <- ccl23.ortho %>% select(-V1)
  colnames(ccl23.ortho) <- c("ccl23_ortho")
  
  ccl24.ortho <- read.table("02_ck_seq/data/processed/ccl24_ortho.txt", header = FALSE)
  ccl24.ortho <- ccl24.ortho %>% select(-V1)
  colnames(ccl24.ortho) <- c("ccl24_ortho")
  
  ccl25.ortho <- read.table("02_ck_seq/data/processed/ccl25_ortho.txt", header = FALSE)
  ccl25.ortho <- ccl25.ortho %>% select(-V1)
  colnames(ccl25.ortho) <- c("ccl25_ortho")
  
  ccl26.ortho <- read.table("02_ck_seq/data/processed/ccl26_ortho.txt", header = FALSE)
  ccl26.ortho <- ccl26.ortho %>% select(-V1)
  colnames(ccl26.ortho) <- c("ccl26_ortho")
  
  ccl27.ortho <- read.table("02_ck_seq/data/processed/ccl27_ortho.txt", header = FALSE)
  ccl27.ortho <- ccl27.ortho %>% select(-V1)
  colnames(ccl27.ortho) <- c("ccl27_ortho")
  
  ccl28.ortho <- read.table("02_ck_seq/data/processed/ccl28_ortho.txt", header = FALSE)
  ccl28.ortho <- ccl28.ortho %>% select(-V1)
  colnames(ccl28.ortho) <- c("ccl28_ortho")
  
  cxcl1.ortho <- read.table("02_ck_seq/data/processed/cxcl1_ortho.txt", header = FALSE)
  cxcl1.ortho <- cxcl1.ortho %>% select(-V1)
  colnames(cxcl1.ortho) <- c("cxcl1_ortho")
  
  cxcl2.ortho <- read.table("02_ck_seq/data/processed/cxcl2_ortho.txt", header = FALSE)
  cxcl2.ortho <- cxcl2.ortho %>% select(-V1)
  colnames(cxcl2.ortho) <- c("cxcl2_ortho")
  
  cxcl3.ortho <- read.table("02_ck_seq/data/processed/cxcl3_ortho.txt", header = FALSE)
  cxcl3.ortho <- cxcl3.ortho %>% select(-V1)
  colnames(cxcl3.ortho) <- c("cxcl3_ortho")
  
  cxcl4.ortho <- read.table("02_ck_seq/data/processed/cxcl4_ortho.txt", header = FALSE)
  cxcl4.ortho <- cxcl4.ortho %>% select(-V1)
  colnames(cxcl4.ortho) <- c("cxcl4_ortho")
  
  cxcl4l1.ortho <- read.table("02_ck_seq/data/processed/cxcl4l1_ortho.txt", header = FALSE)
  cxcl4l1.ortho <- cxcl4l1.ortho %>% select(-V1)
  colnames(cxcl4l1.ortho) <- c("cxcl4l1_ortho")
  
  cxcl5.ortho <- read.table("02_ck_seq/data/processed/cxcl5_ortho.txt", header = FALSE)
  cxcl5.ortho <- cxcl5.ortho %>% select(-V1)
  colnames(cxcl5.ortho) <- c("cxcl5_ortho")
  
  cxcl6.ortho <- read.table("02_ck_seq/data/processed/cxcl6_ortho.txt", header = FALSE)
  cxcl6.ortho <- cxcl6.ortho %>% select(-V1)
  colnames(cxcl6.ortho) <- c("cxcl6_ortho")
  
  cxcl7.ortho <- read.table("02_ck_seq/data/processed/cxcl7_ortho.txt", header = FALSE)
  cxcl7.ortho <- cxcl7.ortho %>% select(-V1)
  colnames(cxcl7.ortho) <- c("cxcl7_ortho")
  
  cxcl8.ortho <- read.table("02_ck_seq/data/processed/cxcl8_ortho.txt", header = FALSE)
  cxcl8.ortho <- cxcl8.ortho %>% select(-V1)
  colnames(cxcl8.ortho) <- c("cxcl8_ortho")
  
  cxcl9.ortho <- read.table("02_ck_seq/data/processed/cxcl9_ortho.txt", header = FALSE)
  cxcl9.ortho <- cxcl9.ortho %>% select(-V1)
  colnames(cxcl9.ortho) <- c("cxcl9_ortho")
  
  cxcl10.ortho <- read.table("02_ck_seq/data/processed/cxcl10_ortho.txt", header = FALSE)
  cxcl10.ortho <- cxcl10.ortho %>% select(-V1)
  colnames(cxcl10.ortho) <- c("cxcl10_ortho")
  
  cxcl11.ortho <- read.table("02_ck_seq/data/processed/cxcl11_ortho.txt", header = FALSE)
  cxcl11.ortho <- cxcl11.ortho %>% select(-V1)
  colnames(cxcl11.ortho) <- c("cxcl11_ortho")
  
  cxcl12.ortho <- read.table("02_ck_seq/data/processed/cxcl12_ortho.txt", header = FALSE)
  cxcl12.ortho <- cxcl12.ortho %>% select(-V1)
  colnames(cxcl12.ortho) <- c("cxcl12_ortho")
  
  cxcl13.ortho <- read.table("02_ck_seq/data/processed/cxcl13_ortho.txt", header = FALSE)
  cxcl13.ortho <- cxcl13.ortho %>% select(-V1)
  colnames(cxcl13.ortho) <- c("cxcl13_ortho")
  
  cxcl14.ortho <- read.table("02_ck_seq/data/processed/cxcl14_ortho.txt", header = FALSE)
  cxcl14.ortho <- cxcl14.ortho %>% select(-V1)
  colnames(cxcl14.ortho) <- c("cxcl14_ortho")
  
  cxcl16.ortho <- read.table("02_ck_seq/data/processed/cxcl16_ortho.txt", header = FALSE)
  cxcl16.ortho <- cxcl16.ortho %>% select(-V1)
  colnames(cxcl16.ortho) <- c("cxcl16_ortho")
  
  cxcl17.ortho <- read.table("02_ck_seq/data/processed/cxcl17_ortho.txt", header = FALSE)
  cxcl17.ortho <- cxcl17.ortho %>% select(-V1)
  colnames(cxcl17.ortho) <- c("cxcl17_ortho")
  
  cx3cl1.ortho <- read.table("02_ck_seq/data/processed/cx3cl1_ortho.txt", header = FALSE)
  cx3cl1.ortho <- cx3cl1.ortho %>% select(-V1)
  colnames(cx3cl1.ortho) <- c("cx3cl1_ortho")
  
  # combine
  data <- bind_cols(all.para, all.cc.cxc.para, cc.para, cxc.para,
                    ccl1.ortho, ccl2.ortho, ccl3.ortho, ccl3l1.ortho, ccl4.ortho, 
                    ccl5.ortho, ccl7.ortho, ccl8.ortho, ccl11.ortho, ccl13.ortho,
                    ccl14.ortho, ccl15.ortho, ccl16.ortho, ccl17.ortho,
                    ccl18.ortho, ccl19.ortho, ccl20.ortho, ccl21.ortho, ccl22.ortho,
                    ccl23.ortho, ccl24.ortho, ccl25.ortho, ccl26.ortho, ccl27.ortho,
                    ccl28.ortho, cxcl1.ortho, cxcl2.ortho, cxcl3.ortho, cxcl4.ortho,
                    cxcl4l1.ortho, cxcl5.ortho, cxcl6.ortho, cxcl7.ortho, cxcl8.ortho,
                    cxcl9.ortho, cxcl10.ortho, cxcl11.ortho, cxcl12.ortho,
                    cxcl13.ortho, cxcl14.ortho, cxcl16.ortho, cxcl17.ortho, cx3cl1.ortho)
  
  rm(all.para, all.cc.cxc.para, cc.para, cxc.para,
     ccl1.ortho, ccl2.ortho, ccl3.ortho, ccl3l1.ortho, ccl4.ortho, 
     ccl5.ortho, ccl7.ortho, ccl8.ortho, ccl11.ortho, ccl13.ortho,
     ccl14.ortho, ccl15.ortho, ccl16.ortho, ccl17.ortho,
     ccl18.ortho, ccl19.ortho, ccl20.ortho, ccl21.ortho, ccl22.ortho,
     ccl23.ortho, ccl24.ortho, ccl25.ortho, ccl26.ortho, ccl27.ortho,
     ccl28.ortho, cxcl1.ortho, cxcl2.ortho, cxcl3.ortho, cxcl4.ortho,
     cxcl4l1.ortho, cxcl5.ortho, cxcl6.ortho, cxcl7.ortho, cxcl8.ortho,
     cxcl9.ortho, cxcl10.ortho, cxcl11.ortho, cxcl12.ortho,
     cxcl13.ortho, cxcl14.ortho, cxcl16.ortho, cxcl17.ortho, cx3cl1.ortho)

  # assign column names
  aln.names <- c(read.table("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
  aln.names <- as.data.frame(aln.names)
  aln.names <- t(aln.names)
  aln.names <- as.data.frame(aln.names)
  colnames(aln.names) <- c("ccn")
  data <- bind_cols(aln.names, data)
  rownames(data) <- 1:nrow(data)
  rm(aln.names)
  
  # spread data
  data <- data %>% gather(protein, ortho_cons, 6:ncol(data))
  
  # validate that scores are accurately assigned to columns
  # N-term checks: checked 0's appropriately assigned where all gaps
  
  # (2.2) ADD ADDITIONAL DESCRIPTORS -------------------------------------------
  # add domain
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_gnccn_to_domain.csv")
  lookup <- lookup %>% unique()
  data$dom <- lookup$dom[match(unlist(data$ccn), lookup$bwccn)]
  rm(lookup)
  
  # add gap numbers (CKR orthologs; 951 seqs) and percentages
  aln <- readAAMultipleAlignment("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.fasta")
  aln <- as.matrix(aln)
  aln.names <- c(read.table("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
  colnames(aln) <- as.matrix(aln.names)
  aln <- as.data.frame(aln)
  #rm(aln.names)
  
  gaps <- NULL
  for (i in 1:ncol(aln)){
    temp <- as.data.frame(table(aln[,i]))
    colnames(temp) <- c("resno", "n")
    temp <- subset(temp, resno == "-")
    temp$pos <- i
    gaps <- rbind.data.frame(gaps, temp)
    rm(temp)
  }
  rm(i)
  gaps$pct <- gaps$n/nrow(aln)
  aln.names <- t(as.data.frame(aln.names))
  colnames(aln.names) <- c("ccn")
  rm(aln)
  
  gaps <- bind_cols(aln.names, gaps)
  rm(aln.names)
  colnames(gaps) <- c("ccn", "restype", "n", "pos", "pct")
  gaps <- select(gaps, ccn, n, pct)
  gaps$ccn <- as.character(gaps$ccn)
  
  data <- left_join(data, gaps)
  colnames(data)[9] <- c("ngaps")
  colnames(data)[10] <- c("pctgaps")
  rm(gaps)
  
  
  # (2.3) REORDER, WRITE OUTPUT ------------------------------------------------
  data <- data %>% select(protein, ccn, dom, all_para, all_cc_cxc_para,
                          cc_para, cxc_para, ortho_cons, ngaps, pctgaps)
  
  # write_csv(data, "02_ck_seq/output/CK_CONSERVATION.csv") LAST WRITTEN 20210323
  