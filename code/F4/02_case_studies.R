# Name:     02_pairwise_rin_compare.R
# Updated:  20230126
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### FUNCTIONS ################################################################

# NOTE REDONE IN FINAL VERSION - HAS NOT CHANGED FROM PRIOR

GetSharedRINs <- function(PDB1, PDB2){
  a <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(file == PDB1) %>%
    select(source_gnccn, target_gnccn)
  b <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(file == PDB2) %>%
    select(source_gnccn, target_gnccn)
  a.b <- intersect(a, b)
  return(a.b)
  rm(a, b, a.b)
}

GetUniqueRINs <- function(PDB1, PDB2){
  a <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(file == PDB1) %>%
    select(source_gnccn, target_gnccn)
  b <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(file == PDB2) %>%
    select(source_gnccn, target_gnccn)
  a.b <- setdiff(a, b)
  return(a.b)
  rm(a, b, a.b)
}


################################################################################

# (1) CXCL12:CXCR4 (ngo) vs. CCL15:CCR1 (7vl9)
shared.ngo.7vl9 <- GetSharedRINs("ngo", "7vl9")
# only.ngo.7vl9 <- GetUniqueRINs("ngo", "7vl9")
# only.7vl9.ngo <- GetUniqueRINs("7vl9", "ngo")
# B3.3-NTr.Cm1 (T1)
# Handful of others that are of the overalapping RIN but unique sequence type

# (2) CCL15:CCR1 (7vl9) vs. CCL20:CCR6 (6wwz)
shared.7vl9.6wwz <- GetSharedRINs("7vl9", "6wwz")
# only.7vl9.6wwz <- GetUniqueRINs("7vl9", "6wwz")
# only.6wwz.7vl9 <- GetUniqueRINs("6wwz", "7vl9")
# 13 overlapping
# CX.1-NTr.Cm1 (T1)
# NTc.Cm3-1x28 (T2)

# (3) CCL3:CCR5 (7f1t) vs CCL15:CCR1 (7vl9)
shared.7f1t.7vl9 <- GetSharedRINs("7f1t", "7vl9")
# only.7f1t.7vl9 <- GetUniqueRINs("7f1t", "7vl9")
# only.7vl9.7f1t <- GetUniqueRINs("7vl9", "7f1t")
# 37 overlapping
# CX.1-NTr.Cm1 (T1)
# CX.1-1x22 (T1)
# NTc.Cm3-1x28 (T2)
# many B1-ECL2
# many b1b2-TM5

# (4) CCL3:CCR5 (7f1t) vs CCL5[5P7]:CCR5 (5uiw)
shared.7f1t.5uiw <- GetSharedRINs("7f1t", "5uiw")
# only.7f1t.5uiw <- GetUniqueRINs("7f1t", "5uiw")
# only.5uiw.7f1t <- GetUniqueRINs("5uiw", "7f1t")
# 53 overlapping
# CX.1-NTr.Cm1 (T1)
# CX.1-1x22 (T1)
# B3.3-NTr.Cm1
# NTc.Cm3-1x28 (T2)
# many B1-ECL2
# many b1b2-TM5
# many B3-NTr

# TESTING - CCL3:CCR5 (7f1t) vs. CCL20:CCR6 (6wwz)
shared.7f1t.6wwz <- GetSharedRINs("7f1t", "6wwz")
only.7f1t.6wwz <- GetUniqueRINs("7f1t", "6wwz")
only.6wwz.7f1t <- GetUniqueRINs("6wwz", "7f1t")
