source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(!(file %in% c("6meo"))) %>%
  filter(no_pdb == 1 ) %>%
  select(source_gnccn, target_gnccn, dom1, dom2) %>%
  unique() 

paste0("There are ", nrow(data), " contacts among 16 chemokine-GPCR complexes",
       " that are unique to a single complex")
