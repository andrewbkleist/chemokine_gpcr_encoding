# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)

################################################################################

data <- read_csv("10_structure_pairwise/output/pairwise_complex_rmsd.csv")
temp <- data %>% count(gnccn) %>% filter(n > 55) 
temp <- temp$gnccn
  # count number of pairwise comparisons in which each position is represented;
  # remove CCN found in less than 50% pairwise comparisons (below)
data <- data %>% group_by(gnccn) %>% 
  filter(gnccn %in% c(temp)) %>%
  mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  select(-RMSD, -file1, -file2) %>% unique() 
rm(temp)

# import pdb, map CCN numbers
pdb <- read.pdb("01_structure_contacts/data/pdbs/7f1t_clean.pdb")
pdb.atom <- pdb$atom
lookup <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20220729.csv")
pdb.atom$source_gnccn <- lookup$ccn_7f1t_ck[match(unlist(pdb.atom$resno), lookup$clean_7f1t_ck)]
rm(lookup)

pdb.atom <- pdb.atom %>%
  mutate(source_gnccn = case_when(
    chain == "A" ~ source_gnccn,
    chain == "B" ~ "NA"
  ))

# map values to pdb b factor
pdb.atom$b <- 0
pdb.atom$b <- data$mean[match(unlist(pdb.atom$source_gnccn), data$gnccn)]
pdb.atom <- pdb.atom %>% mutate(b = case_when(
  is.na(b) ~ mean(data$mean/2),
  !is.na(b) ~ b
))

# write pdb
write.pdb(pdb=pdb, b = pdb.atom$b, file = "F2S/output/5uiw_ca_rmsd.pdb")

