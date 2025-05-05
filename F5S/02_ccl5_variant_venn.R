# Updated:  20201125
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# import contacts
rin.5uiw <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(file == "5uiw") %>%
  select(source_gnccn, target_gnccn) %>% unique()

rin.7o7f <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(file == "7o7f") %>%
  select(source_gnccn, target_gnccn) %>% unique()

rin.5uiw.not.7o7f <- setdiff(rin.5uiw, rin.7o7f)
rin.7o7f.not.5uiw <- setdiff(rin.7o7f, rin.5uiw)
rin.7o7f.and.5uiw <- intersect(rin.7o7f, rin.5uiw)
# # (1) 5uiw not 7o7f
# rin.5uiw.not.7o7f <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
#                               rin.5uiw.not.7o7f,
#                               "5uiw",
#                               "01_structure_contacts/data/pdbs/5uiw_ck_clean.pdb",
#                               "50_rin_alone/output/5uiw_not_7o7f_rins.csv")
# 
# # (2) 7o7f not 5uiw
# rin.7o7f.not.5uiw <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
#                               rin.7o7f.not.5uiw,
#                               "7o7f",
#                               "01_structure_contacts/data/pdbs/7o7f_clean.pdb",
#                               "50_rin_alone/output/7o7f_not_5uiw_rins.csv")
