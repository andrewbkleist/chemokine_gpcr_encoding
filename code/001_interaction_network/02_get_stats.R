# Imports "Supplementary_Table_1.csv" and derives summary information (e.g.
# number of interactions, number of unique sources, etc). Prints to command
# line in R, output stats are referenced in text of manuscript

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
data <- read_csv("data/network/Supplementary_Table_1.csv")
paste0("The number of unique chemokine-GPCR-evidence entries is ", nrow(data))

# pairings considered
ck.gpcr.pairs <- data %>% select(chemokine, gpcr) %>% unique()
paste0("The number of unique chemokine-GPCR PAIRINGS considered (including POSITIVE & NEGATIVE interactions) is ", 
       nrow(ck.gpcr.pairs), " of a total of ", 46*23, 
       " possible pairings, or ",
       round(nrow(ck.gpcr.pairs)/(46*23)*100, digits = 0),
       "% of all possible chemokine-GPCR pairings")
rm(ck.gpcr.pairs)

# positive pairings
pos.ck.gpcr <- data %>% 
  select(chemokine, gpcr, ligand_type, interaction_strength, evidence_grade) %>%
  filter(interaction_strength %in% c(2,3)) %>%
  select(chemokine, gpcr) %>% unique()
paste0("The number of POSITIVE chemokine-GPCR interactions for which evidence was gathered is ",
       nrow(pos.ck.gpcr), " which is ",
       round(nrow(pos.ck.gpcr)/(46*23)*100, digits = 0),
       "% of all possible chemokine-GPCR pairings")
# "positive" interactions are defined as Interaction Strength ≥ 2
# and evidence of any grade
rm(pos.ck.gpcr)

# negative pairings
neg.ck.gpcr <- data %>%
  select(chemokine, gpcr, ligand_type, interaction_strength, evidence_grade) %>%
  filter(interaction_strength %in% c(0))  %>%
  select(chemokine, gpcr) %>% unique()
paste0("The number of NEGATIVE chemokine-GPCR interactions for which evidence was gathered is ",
       nrow(neg.ck.gpcr), " which is ",
       round(nrow(neg.ck.gpcr)/(46*23)*100, digits = 0),
       "% of all possible chemokine-GPCR pairings")
# "positive" interactions are defined as Interaction Strength = 0
# and evidence of any grade
rm(neg.ck.gpcr)

# literature sources
lit <- data %>% select(PMID) %>% unique()
paste0("The number of unique literature SOURCES is ", nrow(lit))
rm(lit)
