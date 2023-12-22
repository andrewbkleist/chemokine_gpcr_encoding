# Imports "Supplementary_Table_1.csv" and selects positive interactions; writes
# output "positive_interactions.csv" file that is opened in Cytoscape
# to generate chemokine network-related figures in manuscript

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
data <- read_csv("data/network/Supplementary_Table_1.csv") %>%
  select(chemokine, gpcr, ligand_type, interaction_strength, evidence_grade) %>%
  filter(interaction_strength %in% c(2,3)) %>%
  select(chemokine, gpcr) %>%
  unique()

# no evidence-documented "high confidence" (i.e. interaction strength ≥ 2)
# pairings for CCL3L1 and CXCL17; will include these with receptor pairings
# as NA in table so that all considered proteins are included in the figure
temp <- data.frame(chemokine = c("CCL3L1", "CXCL17"), gpcr = c(NA, NA))
data <- data %>% rbind(temp)

write_csv(data, "data/network/positive_interactions.csv")


