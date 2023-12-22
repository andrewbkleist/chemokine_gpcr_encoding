source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import contacts
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo") # remove gp120-CCR5

paste0("There are ", 
      nrow(rin), 
      " total contacts among all ", 
      nrow(rin %>% dplyr::select(file) %>% unique()), 
      " complexes")

paste0("There are ", 
       nrow(rin %>% dplyr::select(source_gnccn, target_gnccn) %>% unique()), 
       " unique chemokine-GPCR contacts")
