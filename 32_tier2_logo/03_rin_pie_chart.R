# Name:     03_rin_pie_chart.R
# Updated:  20210330
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
  
################################################################################
  
  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # select and pie chart - CHEMOKINE
  rin %>% #  filter(no_pdb >4) %>% 
    filter(source_gnccn %in% c("b2b3.4")) %>%
    select(source_gnccn, cc_cxc_lr_ck) %>% unique() %>%
    ggplot(aes(x = "", y= cc_cxc_lr_ck)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ source_gnccn)
  
  
  # select and pie chart - RECEPTOR
  rin %>% #  filter(no_pdb >4) %>% 
    filter(target_gnccn %in% c("NTr.Cm1", "1x22", "7x24")) %>%
    select(target_gnccn, all_non_ackr_para_ckr) %>% unique() %>%
    ggplot(aes(x = "", y= all_non_ackr_para_ckr)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    ylim(0,1) +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ target_gnccn)
  
  