source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import data, select only CC/CXC
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0"))

paste0("There are ", nrow(rin), " contacts among CC and CXC chemokine-GPCR complexes")

rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0")) %>%
  filter(cc_cxc_lr_ck >= 0.75, cc_cxc_lr_ckr >= 0.75) %>%
  dplyr::select(source_gnccn, target_gnccn) %>%
  unique()

paste0("There are ", nrow(rin), 
       " distinct contacts among subfamily-predictive residues in CC and CXC chemokine-GPCR complexes")
