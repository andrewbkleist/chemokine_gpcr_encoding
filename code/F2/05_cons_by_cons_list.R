source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import data
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo")

rin <- rin %>% 
  filter((all_cc_cxc_para_ck >= 0.5) & (all_non_ackr_para_ckr >= 0.5))  %>%
  dplyr::select(file, ck, ckr, source_gnccn, target_gnccn, no_pdb)

paste0("The number of files with at least one conserved contact is ",
       nrow(rin %>% dplyr::select(file) %>% unique()))

paste0("The number of coserved contacts is ",
       nrow(rin %>% dplyr::select(source_gnccn, target_gnccn) %>% unique()))

paste0("The coserved contacts are: ")
(rin %>% dplyr::select(source_gnccn, target_gnccn) %>% unique())


