# Name:     02_get_cons_contacts.R
# Updated:  20201115
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1:  GET T2-T2  CHEMOKINE-GPCR CONTACTS ###################################

  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")

  # plot
  rin %>%
    ggplot(aes(cc_cxc_lr_ckr, cc_cxc_lr_ck, size = no_pdb))  +
    geom_point(shape = 21, colour = "black", fill = "grey50",stroke = 0.3) +
    scale_size(range = c(2, 8)) +
    xlim(0.5,1) +
    ylim(0.5,1) +
    theme_minimal()
  rin <- rin %>% filter(all_ortho_ck > 0.4 & all_ortho_ckr > 0.6)
  