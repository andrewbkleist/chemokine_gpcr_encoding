# Name:     07_cxcr6_functional.R
# Updated:  20210712
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

  data <- read_csv("60_snp_cancer/data/processed/cxcr6_functional.csv")
  
  order = c("wt", "e3q")
  order2 = c("pct_wt_binding", "pct_wt_adhesion", "pct_wt_surface_expression")
  data$mutant <- factor(data$mutant, levels = order)
  data$assay <- factor(data$assay, levels = order2)
  
  data %>%
    ggplot(aes(mutant, value))+
    geom_bar(stat = "identity") +
    theme_minimal() +
    facet_grid(. ~ assay)
