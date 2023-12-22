# NOTE: no evidence-documented "high confidence" (i.e. interaction strength ≥ 2)
# pairings for CCL3L1 and CXCL17; these were included with NA receptor pairings
# in "data/network/positive_interactions.csv"

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/network/positive_interactions.csv") %>%
  filter(grepl("CC", chemokine) | grepl("CXC", chemokine)) %>%
  filter(grepl("CC", gpcr) | grepl("CXC", gpcr))

# write_csv(data, "output/F3/cc_cxc_interactions.csv")
# WRITTEN 20231106
# OPENED IN CYTOSCAPE AND NETWORK FOR FIGURE GENERATED THERE
