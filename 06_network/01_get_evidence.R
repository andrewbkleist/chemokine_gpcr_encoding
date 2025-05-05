# Name:     01_get_evidence.R
# Updated:  20201103
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: SELECT CK-CKR INTERACTIONS WITH EVIDENCE >=3 ##########################

  data <- read_csv("06_network/data/raw/ck_ckr_network.csv") %>%
    filter(evidence >=4) %>% select(ck, ckr) %>% unique() %>% 
    filter(ckr != "CXCR8") %>% filter(ckr != "H4 Receptor")
  
  write_csv(data, "06_network/output/chemokine_gpcr_interactome_validated.csv")
  