# Name:     01_cons_by_cons_plot.R
# Updated:  20221020
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE 1 (20210407):
# Note that ECL2 residues prior to 45x50 have poor structural conservation
# due to variability in length of the loop between the two ECL2 beta sheets;
# consequently, changing the alignment to an ECL2.Cm nomenclature by removing
# gaps prior to 45x50 would not have approximated ECL2 structures and would
# thus have little meaning. As  such, ECL2 conservation scores prior to
# 45x50 are not comparable across paralogs and are excluded from the contacts
# represented in the figure

##### 1:  GET CONSERVED  CHEMOKINE-GPCR CONTACTS ###############################

  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # plot
  rin %>%
    # filter(no_pdb >4) %>%
    ggplot(aes(all_non_ackr_para_ckr, all_para_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
    # scale_size(range = c(2, 8)) +
    xlim(0,1) +
    ylim(0,1) +
    theme_minimal()
  
  
  # rin <- rin %>% filter(all_ortho_ck > 0.4 & all_ortho_ckr > 0.6)
  rin <- rin %>% filter((all_cc_cxc_para_ck > 0.5) & (all_non_ackr_para_ckr > 0.5))  %>%
     select(file, ck, ckr, source_gnccn, target_gnccn, no_pdb) %>% unique() %>%
     select(file)  %>% unique()
  
  
  