# Name:     03_vmip_ccl5_pie.R
# Updated:  20230108
# Author:   Andrew Kleist

# packages, working directsory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

  # select PDBs of interest
  pdbs <- c("zheng", "4rws")

  # select tier 2 positions and then select position probability scores
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(cc_cxc_lr_ck >= 0.75) %>%
    filter(file %in% pdbs) %>%
    select(ck, source_gnccn, cc_cxc_lr_score_ck) %>% 
    unique()
  
  # pie chart - add label
  rin <- rin %>% mutate(cc_cxc = case_when(
    cc_cxc_lr_score_ck > 0.5 ~ "cxc",
    cc_cxc_lr_score_ck < 0.5 ~ "cc"
  ))
  rin <- rin %>% count(ck, cc_cxc)
  rin <- rin %>% mutate(pct = case_when(
    ck == "vmipii" ~ n/9,
    ck == "ccl5" ~ n/22
    ))
  
  # pie chart plot
  rin %>%
    ggplot(aes(x = "", y= pct)) +
    geom_bar(width = 1,size = 1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ ck)
  
  