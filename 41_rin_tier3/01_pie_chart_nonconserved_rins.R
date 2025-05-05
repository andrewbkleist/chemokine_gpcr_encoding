# Name:     01_cons_by_cons_plot.R
# Updated:  20210406
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE 1 (20210406):
# Note that ECL2 residues prior to 45x50 have poor structural conservation
# due to variability in length of the loop between the two ECL2 beta sheets;
# consequently, changing the alignment to an ECL2.Cm nomenclature by removing
# gaps prior to 45x50 would not have approximated ECL2 structures and would
# thus have little meaning. As  such, ECL2 conservation scores prior to
# 45x50 are not comparable across paralogs and are excluded from the contacts
# represented in the figure

##### 1:  GET NON-CONSERVED  CHEMOKINE-GPCR CONTACTS ############################

  # import data
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # classify based on quadrant
  rin <- rin %>% mutate(all_para_quad = case_when(
    all_para_ck < 0.5 & all_non_ackr_para_ckr < 0.5 ~ "BL",
    all_para_ck < 0.5 & all_non_ackr_para_ckr > 0.5 ~ "BR",
    all_para_ck > 0.5 & all_non_ackr_para_ckr < 0.5 ~ "TL",
    all_para_ck > 0.5 & all_non_ackr_para_ckr > 0.5 ~ "TR"
  )) %>%
    filter(!is.na(all_para_quad))
  
  # count number per quardant
  no <- rin %>% count(all_para_quad)
  no <- no %>% mutate(pct_all_para_quad = n / nrow(rin))
  
  
  # PIE CHART ALL
  no %>%
    ggplot(aes(x = "", y= pct_all_para_quad)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # PIE CHART EXAMPLES ALL FOUR...
  
  
  
  # rin <- rin %>% filter(all_ortho_ck < 0.5 & all_ortho_ckr < 0.5)
  test <- rin %>% filter((all_para_ck < 0.5) & (all_non_ackr_para_ckr < 0.5))  %>%
    filter((ortho_cons_ck > 0.5) & (ortho_cons_ckr > 0.5)) # %>%
    #select(source_gnccn, target_gnccn, no_pdb) %>% unique()
  
  
  # plot (with sizing by no_pdb)
  rin %>%
       # filter(no_pdb >4) %>%
       ggplot(aes(all_non_ackr_para_ckr, all_para_ck, size = no_pdb))  +
       geom_point(shape = 21, colour = "black", fill = "grey70",stroke = 0.3) +
       scale_size(range = c(2, 8)) +
       xlim(0,1) +
       ylim(0,1) +
       theme_minimal()