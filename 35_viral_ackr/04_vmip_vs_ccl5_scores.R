# Name:     04_vmip_vs_ccl5_scores.R
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
    select(ck, source_gnccn, cc_cxc_lr_score_ck, cc_cxc_lr_score_sd_ck) %>% 
    unique()
  
  vmipii.int <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(cc_cxc_lr_ck >= 0.75) %>%
    filter(file == "4rws") %>%
    select(source_gnccn) %>% 
    unique()
  
  vmipii.int <- vmipii.int$source_gnccn
  
  # # select only residues for which vMIP-II has contacts
  # rin <- rin %>% select(-cc_cxc_lr_score_sd_ck) %>% pivot_wider(names_from = ck, values_from = cc_cxc_lr_score_ck)
  # rin <- rin %>% filter(!(is.na(vmipii))) %>% filter(!(is.na(ccl5)))
  # rin <- rin %>% mutate(diff = vmipii-ccl5)
  rin <- rin %>% filter(source_gnccn %in% vmipii.int)

  # order
  order.ccn <- as.factor(c("NTc.Cm10","NTc.Cm9","NTc.Cm8","NTc.Cm7","NTc.Cm6","NTc.Cm5","NTc.Cm4","NTc.Cm3","NTc.Cm2","NTc.Cm1", "CX.1","CX.2","CX.3","CX.4","CX.5","cxb1.1","cxb1.2","cxb1.3","cxb1.4","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.9","cxb1.10","cxb1.11","cxb1.12","cxb1.13","cxb1.14","cxb1.15","cxb1.16","cxb1.17","cxb1.18","cxb1.19","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.1","b1b2.2","b1b2.3","b1b2.4","b1b2.5","b1b2.6","b1b2.7","b1b2.8","b1b2.9","b1b2.10","b1b2.11","b1b2.12","b1b2.13","b1b2.14","b1b2.15","b1b2.16","b1b2.17","b1b2.18","b1b2.19","b1b2.20","b1b2.21","b1b2.22","b1b2.23","b1b2.24","b1b2.25","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.1","b2b3.2","b2b3.3","b2b3.4","b2b3.5","b2b3.6","b2b3.7","b2b3.8","b2b3.9","b2b3.10","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","b3h.5","b3h.6","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10"))
  rin$source_gnccn <- factor(rin$source_gnccn, levels = order.ccn)
  
  # plot
  rin %>%
    filter(source_gnccn %in% c("NTc.Cm4", "NTc.Cm1", "cxb1.1")) %>%
    ggplot(aes(source_gnccn, cc_cxc_lr_score_ck, fill = ck)) +
    geom_errorbar(aes(ymin=cc_cxc_lr_score_ck-cc_cxc_lr_score_sd_ck, ymax=cc_cxc_lr_score_ck+cc_cxc_lr_score_sd_ck), width=.2) +
    geom_point(shape = 21,   size = 4, stroke = 0.5) +
    theme_minimal() +
    # ylim(0,1) +
    scale_fill_manual(values=c("steelblue4", "red4")) +
    coord_flip() +
    theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0))
  
  