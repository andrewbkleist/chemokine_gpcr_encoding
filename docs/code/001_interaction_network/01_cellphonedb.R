# Imports .csv files containing raw ligand-receptor pairs from CellphoneDB 
# (https://www.cellphonedb.org) and merges into single file, removes redundant
# interactions, and removes non chemokine/chemokine GPCR-containing interactions
# source libraries, functions. Sources removed, note that majority reference PMID
# 24218476 (IUPHAR chemokine/GPCR guide paper); only CXCL14-CXCR4 (PMID 32212206)
# and CCL3-CCR5 (PMID 15251452) reference specific, primary research articles.
# Data downloaded from Cellphonedb on 08.11.2023.

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
a <- read_csv("data/network/cellphonedb/raw/ligand_set_1.csv")
b <- read_csv("data/network/cellphonedb/raw/ligand_set_2.csv")
c <- read_csv("data/network/cellphonedb/raw/ligand_set_3.csv")
d <- read_csv("data/network/cellphonedb/raw/ligand_set_4.csv")
e <- read_csv("data/network/cellphonedb/raw/gpcr_set_1.csv")
f <- read_csv("data/network/cellphonedb/raw/gpcr_set_2.csv")

data <- do.call("rbind", list(a,b,c,d,e,f))
rm(a,b,c,d,e,f)
data <- data %>%
  filter(! (`Gene name A` %in% c("RARRES2"))) %>%
  mutate(`Gene name A` = case_when(
    `Gene name A` == "PF4" ~ "CXCL4",
    `Gene name A` == "PPBP" ~ "CXCL7",
    !(`Gene name A` %in% c("PF4", "PPBP")) ~ `Gene name A`
  )) %>%
  filter(! (`Gene name B` %in% c("DPP4","HRH4","FPR2"))) %>%
  unique() %>%
  select(`Gene name A`,`Gene name B`)
colnames(data)[1] <- c("chemokine")
colnames(data)[2] <- c("gpcr")
write_csv(data, "data/network/cellphonedb/clean/cellphonedb_chemokine_gpcr.csv")
  


