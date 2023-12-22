source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/network/Supplementary_Table_1.csv") %>%
  dplyr::select(chemokine, gpcr, author_date_journal, PMID, binding_kd_ki_nm, 
         binding_ec50_ic50_nm, signaling_ec50_ic50_nm, chemotaxis_ec50_ic50_nm, 
         chemotaxis_max_nm, ligand_type, interaction_strength, evidence_grade)

# unique entries
cat(paste0("- ", "#### The number of unique chemokine-GPCR-evidence entries is **", nrow(data), "**"))

# pairings considered
ck.gpcr.pairs <- data %>% dplyr::select(chemokine, gpcr) %>% unique()
cat(paste0("- ", "#### The number of unique chemokine-GPCR **pairings** considered (including **positive** & **negative** interactions) is ", 
           nrow(ck.gpcr.pairs), " of a total of ", 46*23, 
           " possible pairings, or **",
           round(nrow(ck.gpcr.pairs)/(46*23)*100, digits = 0),
           "%** of all possible chemokine-GPCR pairings"))
rm(ck.gpcr.pairs)


# positive pairings
pos.ck.gpcr <- data %>% 
  dplyr::select(chemokine, gpcr, ligand_type, interaction_strength, evidence_grade) %>%
  filter(interaction_strength %in% c(2,3)) %>%
  dplyr::select(chemokine, gpcr) %>% unique()
cat(paste0("- ", "#### The number of **positive** chemokine-GPCR interactions for which evidence was gathered is ",
           nrow(pos.ck.gpcr), " which is **",
           round(nrow(pos.ck.gpcr)/(46*23)*100, digits = 0),
           "%** of all possible chemokine-GPCR pairings"))
# "positive" interactions are defined as Interaction Strength ≥ 2
# and evidence of any grade
rm(pos.ck.gpcr)

# negative pairings
neg.ck.gpcr <- data %>%
  dplyr::select(chemokine, gpcr, ligand_type, interaction_strength, evidence_grade) %>%
  filter(interaction_strength %in% c(0))  %>%
  dplyr::select(chemokine, gpcr) %>% unique()
cat(paste0("- ", "#### The number of **negative** chemokine-GPCR interactions for which evidence was gathered is ",
           nrow(neg.ck.gpcr), " which is **",
           round(nrow(neg.ck.gpcr)/(46*23)*100, digits = 0),
           "%** of all possible chemokine-GPCR pairings"))
# "positive" interactions are defined as Interaction Strength = 0
# and evidence of any grade
rm(neg.ck.gpcr)

# literature sources
lit <- data %>% dplyr::select(PMID) %>% unique()
cat(paste0("- ", "#### The number of unique literature **sources** is **", nrow(lit), "**"))
rm(lit)
