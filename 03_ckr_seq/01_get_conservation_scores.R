# Name:     01_get_conservation_scores.R
# Updated:  20210328
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

# NOTE 1 (20200919):
# TOPIC: ECL2
# Aligning CKRs, one can make argument that ECL2b can be considered structuarlly 
# conserved (see /data/pdb/ecl2.pse) and thus alignable
# Note that current alignment from which all calculations are made is called
# FULL_RECEPTOR_ALIGNMENT_CYS_ADJ.fasta which contains all orthologs/paralogs; 
# this file was modified by opening aligment and converting to a data frame,
# writing the ECL2 portion to a sequence - "collapsing" the gaps around the Cys;
# Note that since insertions relative to human have been removed, all 
# modifications were made across the entire set of orthologs column-wise
# for each receptor; ECL2a seems better aligned to TM4 when looking at structures
# (in oart because loop between ECL2 beta sheets is so variable) so will leave alone
# Summary: adjusted region after 4x65 and before 5x31; within this, did not touch
# between ECL2.1 through 45x52, then removed gaps (only shared gaps among receptor-
# specific orthologs) after 45x52 
# Alignment saved as SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta
#
#
# NOTE 2 (20200919):
# TOPIC: NOMENCLATURE
# Updated the nomenclature to reflect new ECL2 nomenclature (note aln length 
# did not change) and changed N-term nomeclature to be exlusively with "NTr.Cp"
# convention; note that ECL2a convention (ECL2.1...) has not changed, only ECL2b,
# since ECL2a is not adjusted relative to Cys
# New nomenclature saved as SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt
#
#
# NOTE 3 (20200919):
# TOPIC: CLASS A ALIGNMENT
# Note that since the ECL2 length has not changed there will be no change
# in the Class A alignment outside of the ECL2 region
# Input fasta: SEQUENCES_CLASSA_AND_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta
# Separated Nterm and ECL2
# Did following...
# (1) NTerm
# - Selected gnNT.Cm50:gn1x30
# - For the following receptors which **have Cys at 7x24** and thus likely to
#   have conserved 7x24-1x22 disulfide, adjusted N-terminal Cys to align to 
#   1x22 position (51st residue position in alignment)
#   agrt1, ap1, bkbr1, bkbr2, ednra, ednrb, cltr1, lpar4, lpar6, gpr55, p2y1,
#   p2y2, p2y4, p2y6, hcar2, hcar3, oxgr1, gpr4, gpr17, gpr25, gpr34, gpr132,
#   gpr171, gpr174, gpr182, gpr183, p2y10
# - For all of above moved 1x23-1x30 positions to approximate 1x22 (ie removed
#   intervening gaps); 1x30 is unaffected in all seqs
# - For the following which have Cys near 7x24 (likely poor alignment), also
#   adjusted Nterm Cys to 1x22 position
#   cltr2, hcar1, sucr1
# (2) ECL2
# - Selected gnECL2.1:gnECL2.Cp11
# - Removed gaps after Cys so all seqs left adjusted to Cys
# - Exceptions include GPCRs for which 4x50 is not Cys (eg lpar, s1pr) which have
#   different ECL2 architecture
# 
# FINAL UPDATED SEQUENCES (all CKR seqs are same in both):
# - SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta             # ckr only
# - SEQUENCES_CLASSA_AND_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta  # class A & ckr
# - SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt               # col names
#
#
# NOTE 4 (20200920):
# TOPIC: ALIGNMENT SETS
# Want the following per-position descriptors:
# - Ortholog conservation (x23)
# - CC paralog conservation (1)
# - CXC paralog conservation (1)
# - CC and CXC paralog conservation (1)
# - ACK paralog conservation (1)
# - All 23 paralog conservation (1)
# - All 951 ortholog conservation (1) - REMOVED, 20210303 NOTE
# - T2 scores (1) [different script]
# - Receptor-specific T2 scores (all class A, x ~250) [different script]
# - Receptor-specific T2 scores (viral CKRs, x 11) [different script]
# **For alignment sets - manually created the different sequence sets
#
# NOTE 5 (20201020):
# TOPIC: VIRAL RECEPTORS
# (1) Sequences from Alex Gunnarsson; includes 11 viral receptor seqs
# (2) Combined 11 viral seqs with 23 human receptor paralogs and aligned using
#     ClustalO: viral_human_ckrs_clustal.fasta
# (3) GPCRdb has US28 sequence - compared Clustal alignment to GPCRdb alignment
#     between US28 and CXCR4; alignments match in core regions
# (4) Pasted 11 viral sequences + CXCR4 into new alignment with 23 human paralogs
#     from the final master alignment; used CXCR4 as "bridge" to adjust viral
#     sequences to match master alignment; ***insertions in viral sequences
#     relative to the master alignment were deleted***; gaps in C-terminus after 
#     NPxxY were remoevd and resiudes left aligned to TM7; ECL2 adjusted as with
#     class A
# (5) Final alignment called SEQUENCES_VIRAL_NTERM_ECL2_UPDATED.fasta and placed
#     in /03_ckr_seq/data/raw; note that it will be considered separately from 
#     class A master alignment for Tier 2 scoring so as not to redo entire T2
#     scoring for these 11 additional sequences
# 
# 
# NOTE 6 (20210303):
# TOPIC: ERROR WITH MSTATX CALC > 500 SEQ
# - Noted that MstatX has 500 seq maximum
# - No other calculations should be affected (ortholog conservation score is
#   the only one using >500 seq)
#
#
##### PART 1: CONSERVATION SCORING #############################################
# From all alignments, calculate per-position conservation scores
# Note that all sequence sets were manually curated
# the following commands were run in the command line from:
# /Users/akleist/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/03_ckr_seq/data/processed
#
# mstatx -i ALL_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_para.txt
# mstatx -i ALL_non_ackr_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_NON_ACKR_para.txt
# mstatx -i ALL_CC_CXC_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_CC_CXC_para.txt
# mstatx -i CC_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o CC_para.txt
# mstatx -i CXC_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o CXC_para.txt
# mstatx -i ACK_para.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ACK_para.txt
# mstatx -i ALL_classa.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ALL_classa.txt
# mstatx -i ccr1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr1_ortho.txt
# mstatx -i ccr2_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr2_ortho.txt
# mstatx -i ccr3_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr3_ortho.txt
# mstatx -i ccr4_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr4_ortho.txt
# mstatx -i ccr5_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr5_ortho.txt
# mstatx -i ccr6_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr6_ortho.txt
# mstatx -i ccr7_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr7_ortho.txt
# mstatx -i ccr8_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr8_ortho.txt
# mstatx -i ccr9_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr9_ortho.txt
# mstatx -i ccr10_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccr10_ortho.txt
# mstatx -i cxcr1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcr1_ortho.txt
# mstatx -i cxcr2_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcr2_ortho.txt
# mstatx -i cxcr3_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcr3_ortho.txt
# mstatx -i cxcr4_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcr4_ortho.txt
# mstatx -i cxcr5_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcr5_ortho.txt
# mstatx -i cxcr6_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cxcr6_ortho.txt
# mstatx -i cx3cr1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o cx3cr1_ortho.txt
# mstatx -i xcr1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o xcr1_ortho.txt
# mstatx -i ackr1_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ackr1_ortho.txt
# mstatx -i ackr2_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ackr2_ortho.txt
# mstatx -i ackr3_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ackr3_ortho.txt
# mstatx -i ackr4_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ackr4_ortho.txt
# mstatx -i ccrl2_ortho.fasta -s trident -m /Applications/MstatX-master/data/aaindex/HENS920102.mat -o ccrl2_ortho.txt

##### PART 2: IMPORT AND TIDY SCORES ###########################################

  # (2.1) IMPORT AND TIDY ------------------------------------------------------
  # import scores
  all.para <- read.table("03_ckr_seq/data/processed/ALL_para.txt", header = FALSE)
  all.para <- all.para %>% select(-V1)
  colnames(all.para) <- c("all_para")
  
  all.non.ackr.para <- read.table("03_ckr_seq/data/processed/ALL_non_ackr_para.txt", header = FALSE)
  all.non.ackr.para <- all.non.ackr.para %>% select(-V1)
  colnames(all.non.ackr.para) <- c("all_non_ackr_para")
  
  all.cc.cxc.para <- read.table("03_ckr_seq/data/processed/ALL_CC_CXC_para.txt", header = FALSE)
  all.cc.cxc.para <- all.cc.cxc.para %>% select(-V1)
  colnames(all.cc.cxc.para) <- c("all_cc_cxc_para")
  
  cc.para <- read.table("03_ckr_seq/data/processed/CC_para.txt", header = FALSE)
  cc.para <- cc.para %>% select(-V1)
  colnames(cc.para) <- c("cc_para")
  
  cxc.para <- read.table("03_ckr_seq/data/processed/CXC_para.txt", header = FALSE)
  cxc.para <- cxc.para %>% select(-V1)
  colnames(cxc.para) <- c("cxc_para")
  
  ack.para <- read.table("03_ckr_seq/data/processed/ACK_para.txt", header = FALSE)
  ack.para <- ack.para %>% select(-V1)
  colnames(ack.para) <- c("ack_para")
  
  # REMOVED
  # all.ortho <- read.table("03_ckr_seq/data/processed/ALL_ortho.txt", header = FALSE)
  # all.ortho <- all.ortho %>% select(-V1)
  # colnames(all.ortho) <- c("all_ortho")
  
  # REMOVED
  # all.cc.cxc.ortho <- read.table("03_ckr_seq/data/processed/ALL_CC_CXC_ortho.txt", header = FALSE)
  # all.cc.cxc.ortho <- all.cc.cxc.ortho %>% select(-V1)
  # colnames(all.cc.cxc.ortho) <- c("all_cc_cxc_ortho")
  
  all.classa <- read.table("03_ckr_seq/data/processed/ALL_classa.txt", header = FALSE)
  all.classa <- all.classa %>% select(-V1)
  colnames(all.classa) <- c("all_classa")
  
  ccr1.ortho <- read.table("03_ckr_seq/data/processed/ccr1_ortho.txt", header = FALSE)
  ccr1.ortho <- ccr1.ortho %>% select(-V1)
  colnames(ccr1.ortho) <- c("ccr1")
  
  ccr2.ortho <- read.table("03_ckr_seq/data/processed/ccr2_ortho.txt", header = FALSE)
  ccr2.ortho <- ccr2.ortho %>% select(-V1)
  colnames(ccr2.ortho) <- c("ccr2")
  
  ccr3.ortho <- read.table("03_ckr_seq/data/processed/ccr3_ortho.txt", header = FALSE)
  ccr3.ortho <- ccr3.ortho %>% select(-V1)
  colnames(ccr3.ortho) <- c("ccr3")
  
  ccr4.ortho <- read.table("03_ckr_seq/data/processed/ccr4_ortho.txt", header = FALSE)
  ccr4.ortho <- ccr4.ortho %>% select(-V1)
  colnames(ccr4.ortho) <- c("ccr4")
  
  ccr5.ortho <- read.table("03_ckr_seq/data/processed/ccr5_ortho.txt", header = FALSE)
  ccr5.ortho <- ccr5.ortho %>% select(-V1)
  colnames(ccr5.ortho) <- c("ccr5")
  
  ccr6.ortho <- read.table("03_ckr_seq/data/processed/ccr6_ortho.txt", header = FALSE)
  ccr6.ortho <- ccr6.ortho %>% select(-V1)
  colnames(ccr6.ortho) <- c("ccr6")
  
  ccr7.ortho <- read.table("03_ckr_seq/data/processed/ccr7_ortho.txt", header = FALSE)
  ccr7.ortho <- ccr7.ortho %>% select(-V1)
  colnames(ccr7.ortho) <- c("ccr7")
  
  ccr8.ortho <- read.table("03_ckr_seq/data/processed/ccr8_ortho.txt", header = FALSE)
  ccr8.ortho <- ccr8.ortho %>% select(-V1)
  colnames(ccr8.ortho) <- c("ccr8")
  
  ccr9.ortho <- read.table("03_ckr_seq/data/processed/ccr9_ortho.txt", header = FALSE)
  ccr9.ortho <- ccr9.ortho %>% select(-V1)
  colnames(ccr9.ortho) <- c("ccr9")
  
  ccr10.ortho <- read.table("03_ckr_seq/data/processed/ccr10_ortho.txt", header = FALSE)
  ccr10.ortho <- ccr10.ortho %>% select(-V1)
  colnames(ccr10.ortho) <- c("ccr10")
  
  cxcr1.ortho <- read.table("03_ckr_seq/data/processed/cxcr1_ortho.txt", header = FALSE)
  cxcr1.ortho <- cxcr1.ortho %>% select(-V1)
  colnames(cxcr1.ortho) <- c("cxcr1")
  
  cxcr2.ortho <- read.table("03_ckr_seq/data/processed/cxcr2_ortho.txt", header = FALSE)
  cxcr2.ortho <- cxcr2.ortho %>% select(-V1)
  colnames(cxcr2.ortho) <- c("cxcr2")
  
  cxcr3.ortho <- read.table("03_ckr_seq/data/processed/cxcr3_ortho.txt", header = FALSE)
  cxcr3.ortho <- cxcr3.ortho %>% select(-V1)
  colnames(cxcr3.ortho) <- c("cxcr3")
  
  cxcr4.ortho <- read.table("03_ckr_seq/data/processed/cxcr4_ortho.txt", header = FALSE)
  cxcr4.ortho <- cxcr4.ortho %>% select(-V1)
  colnames(cxcr4.ortho) <- c("cxcr4")
  
  cxcr5.ortho <- read.table("03_ckr_seq/data/processed/cxcr5_ortho.txt", header = FALSE)
  cxcr5.ortho <- cxcr5.ortho %>% select(-V1)
  colnames(cxcr5.ortho) <- c("cxcr5")
  
  cxcr6.ortho <- read.table("03_ckr_seq/data/processed/cxcr6_ortho.txt", header = FALSE)
  cxcr6.ortho <- cxcr6.ortho %>% select(-V1)
  colnames(cxcr6.ortho) <- c("cxcr6")
  
  cx3cr1.ortho <- read.table("03_ckr_seq/data/processed/cx3cr1_ortho.txt", header = FALSE)
  cx3cr1.ortho <- cx3cr1.ortho %>% select(-V1)
  colnames(cx3cr1.ortho) <- c("cx3cr1")
  
  xcr1.ortho <- read.table("03_ckr_seq/data/processed/xcr1_ortho.txt", header = FALSE)
  xcr1.ortho <- xcr1.ortho %>% select(-V1)
  colnames(xcr1.ortho) <- c("xcr1")
  
  ackr1.ortho <- read.table("03_ckr_seq/data/processed/ackr1_ortho.txt", header = FALSE)
  ackr1.ortho <- ackr1.ortho %>% select(-V1)
  colnames(ackr1.ortho) <- c("ackr1")
  
  ackr2.ortho <- read.table("03_ckr_seq/data/processed/ackr2_ortho.txt", header = FALSE)
  ackr2.ortho <- ackr2.ortho %>% select(-V1)
  colnames(ackr2.ortho) <- c("ackr2")
  
  ackr3.ortho <- read.table("03_ckr_seq/data/processed/ackr3_ortho.txt", header = FALSE)
  ackr3.ortho <- ackr3.ortho %>% select(-V1)
  colnames(ackr3.ortho) <- c("ackr3")
  
  ackr4.ortho <- read.table("03_ckr_seq/data/processed/ackr4_ortho.txt", header = FALSE)
  ackr4.ortho <- ackr4.ortho %>% select(-V1)
  colnames(ackr4.ortho) <- c("ackr4")
  
  ccrl2.ortho <- read.table("03_ckr_seq/data/processed/ccrl2_ortho.txt", header = FALSE)
  ccrl2.ortho <- ccrl2.ortho %>% select(-V1)
  colnames(ccrl2.ortho) <- c("ccrl2")
  
  # combine
  data <- bind_cols(all.para, all.non.ackr.para, all.cc.cxc.para, all.classa, cc.para, cxc.para, ack.para,
                    ccr1.ortho, ccr2.ortho, ccr3.ortho, ccr4.ortho, ccr5.ortho, 
                    ccr6.ortho, ccr7.ortho, ccr8.ortho, ccr9.ortho, ccr10.ortho,
                    cxcr1.ortho, cxcr2.ortho, cxcr3.ortho, cxcr4.ortho,
                    cxcr5.ortho, cxcr6.ortho,
                    ackr1.ortho, ackr2.ortho, ackr3.ortho, ackr4.ortho, ccrl2.ortho,
                    cx3cr1.ortho, xcr1.ortho)
  
  rm(all.para, all.non.ackr.para, all.cc.cxc.para, all.classa, cc.para, cxc.para, ack.para,
     ccr1.ortho, ccr2.ortho, ccr3.ortho, ccr4.ortho, ccr5.ortho, 
     ccr6.ortho, ccr7.ortho, ccr8.ortho, ccr9.ortho, ccr10.ortho,
     cxcr1.ortho, cxcr2.ortho, cxcr3.ortho, cxcr4.ortho,
     cxcr5.ortho, cxcr6.ortho,
     ackr1.ortho, ackr2.ortho, ackr3.ortho, ackr4.ortho, ccrl2.ortho,
     cx3cr1.ortho, xcr1.ortho)
  
  # assign column names
  aln.names <- c(read.table("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ",", colClasses = "character"))
  aln.names <- as.data.frame(aln.names)
  aln.names <- t(aln.names)
  aln.names <- as.data.frame(aln.names)
  colnames(aln.names) <- c("gn")
  data <- bind_cols(aln.names, data)
  rownames(data) <- 1:nrow(data)
  rm(aln.names)
  
  # spread data
  data <- data %>% gather(protein, ortho_cons, 9:ncol(data))
  
  # validate that scores are accurately assigned to columns
  # (1) N-/C- term checks: check  0's appropriately assigned where all gaps:
  #   - CC, CXC, ACK, all_classa, all_para; all_ortho hard to assess
  #   - CCR1-10, CXCR1-6, ACKR1-CCRL2, XCR1, CX3CR1
  
  
  # (2.2) ADD ADDITIONAL DESCRIPTORS -------------------------------------------
  # add domain
  lookup <- read_csv("01_structure_contacts/data/processed/lookup_gnccn_to_domain.csv")
  lookup$gn <- c("gn")
  lookup <- lookup %>% unique()
  lookup <- lookup %>% unite(gn, c(gn, bwccn), sep = "")
  data$dom <- lookup$dom[match(unlist(data$gn), lookup$gn)]
  rm(lookup)
  
  # add gap numbers (CKR orthologs; 951 seqs) and percentages
  aln <- readAAMultipleAlignment("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta")
  aln <- as.matrix(aln)
  aln.names <- c(read.table("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ",", colClasses = "character"))
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
  colnames(aln.names) <- c("gn")
  rm(aln)
  
  gaps <- bind_cols(aln.names, gaps)
  rm(aln.names)
  colnames(gaps) <- c("gn", "restype", "n", "pos", "pct")
  gaps <- select(gaps, gn, n, pct)
  gaps$gn <- as.character(gaps$gn)

  data <- left_join(data, gaps)
  colnames(data)[12] <- c("ngaps")
  colnames(data)[13] <- c("pctgaps")
  rm(gaps)
  
  # (2.3) REORDER, WRITE OUTPUT ------------------------------------------------
  data <- data %>% select(protein, gn, dom, all_para, all_non_ackr_para, all_cc_cxc_para, all_classa, 
                          cc_para, cxc_para, ack_para, ortho_cons, ngaps, pctgaps)
  
  # write_csv(data, "03_ckr_seq/output/CKR_CONSERVATION.csv") LAST WRITTEN 20210328
  
        
  