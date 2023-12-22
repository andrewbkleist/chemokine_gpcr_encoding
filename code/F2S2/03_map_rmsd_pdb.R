source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/structure/processed/pairwise_complex_rmsd.csv")
temp <- data %>% dplyr::count(gnccn) %>% filter(n > 60) 
temp <- temp$gnccn
  # count number of pairwise comparisons in which each position is represented;
  # remove CCN found in less than 50% pairwise comparisons (below)
data <- data %>% group_by(gnccn) %>% 
  filter(gnccn %in% c(temp)) %>%
  dplyr::mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  dplyr::select(-RMSD, -file1, -file2) %>% unique() 
rm(temp)

# import pdb, map CCN numbers
pdb <- read.pdb("data/structure/raw/pdb_files/7f1t_clean.pdb")
pdb.atom <- pdb$atom
lookup <- read_csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv")
pdb.atom$source_gnccn <- lookup$ccn_7f1t_ck[match(unlist(pdb.atom$resno), lookup$clean_7f1t_ck)]
rm(lookup)

pdb.atom <- pdb.atom %>%
  dplyr::mutate(source_gnccn = case_when(
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
# write.pdb(pdb=pdb, b = pdb.atom$b, file = "output/F2S2/5uiw_ca_rmsd.pdb")

# PYLMOL COMMANDS (RUN IN PYMOL) - RANGE TAKES RANGE OF B-FACTORS FROM "data":
# cartoon putty
# spectrum b, firebrick_grey70, minimum=4, maximum=14
