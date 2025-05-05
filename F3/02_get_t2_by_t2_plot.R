# Name:     02_t2_by_t2_plot.R
# Updated:  20221107
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################
  
  # import data, select only CC/CXC
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", "7xa3", "7f1t", "6wwz", "6lfo", "ngo"))
  
  # plot (with sizing by no_pdb)
  rin %>%
    # filter(no_pdb >4) %>%
    ggplot(aes(cc_cxc_lr_ckr, cc_cxc_lr_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
    # scale_size(range = c(2, 8)) +
    xlim(0.5,1) +
    ylim(0.5,1) +
    theme_minimal()
  
  # subset for analysis
  test <- rin %>% filter(cc_cxc_lr_ck >= 0.75 & cc_cxc_lr_ckr >= 0.75) %>%
    count(source_gnccn, target_gnccn) %>% unique()
  
  