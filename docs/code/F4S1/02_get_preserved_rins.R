source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo") %>%
  filter(all_para_ck < 0.5 | all_non_ackr_para_ckr <0.5) %>%
  dplyr::count(source_gnccn, target_gnccn) %>%
  filter(n > 8)
paste0("Preserved contacts (>8/16 complexes) comprised of at least one",
       " nonconserved residue ", "are as follows:")
rin
# note that ECL3.Cm2 is "disqualified" from being a preserved contact since
# ECL2 residues preceding 45x50 are not structurally aligned
  