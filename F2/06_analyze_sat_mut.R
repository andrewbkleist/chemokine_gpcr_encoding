# Name:     01_analyze_sat_mut.R
# Updated:  20221106
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE (20210330):
# To identify mutation sites with known minimal impact on CXCL12 activation
# (ie negative control), we used Wescott, Supp Table 3-5 and identified CXCR4 
# mutations with minimaleffects on signaling as those with ≥95% efficacy; only 
# resi with GPCRdb positions were selected; the list includes...
# 3.53, 3.54, 3.57, 5.64, 6.32, 8.47, 8.49 (Supp Table 3; Gpro interface)
# 5.36(x37), 5.37(x38), 5.40(x41), 5.44(x45) (Supp Table 4; dimerization interface)
# 3.51, 5.59, 5.62, 5.63, 5.64, 6.39, 6.43, 7.49, 7.52 (Supp Table 5; signaling motifs)
# Of these residues, the TM5 ones are at the interface region but do not
# make structural contacts (contacts are 5x33, 5x36, 5x39, 5x40)
#
# As "positive control", from Stephens 2020 Sci Sig paper, 6x58 mutations have 
# most profound effects on CXCL12 binding and are at the interface


##### 1: MUTATIONAL ANALYSIS ###################################################

  # import
  data <- read_csv("09_sat_mut_data/output/cxcr4_clean_gpcrdb_means_log.csv")
  data2 <- data
  
  # plot - BOXPLOT
  data$gn <- factor(data$gn, levels = c("5x37", "5x38", "5x41", "5x45", "6x58", "7x27", "NTr.Cm1", "1x22", "7x24"))
  data %>% filter(gn %in% c("5x38", "6x58",  "NTr.Cm1", "1x22", "7x24")) %>%
    #filter(sub != "\\*") %>%
    select(resno, gn, sele, sub_mean) %>%
    filter(sele == "cxcl12") %>%
    unique() %>%
    ggplot(aes(gn, sub_mean)) +
    geom_boxplot() +
    coord_flip() +
    theme_minimal()
  
  # stats
  test <- data2 %>% filter(gn %in% c("5x38", "NTr.Cm1")) %>% 
    filter(sele == "cxcl12") %>%
    select(resno, gn, value) %>% unique()
  stat.test <- t.test(value ~ gn, data = test, exact = FALSE)
  stat.test$p.value
  
  # stats
  test <- data %>% filter(gn %in% c("5x38", "1x22")) %>% 
    filter(sele == "cxcl12") %>%
    select(resno, gn, value) %>% unique()
  stat.test <- t.test(value ~ gn, data = test, exact = FALSE)
  stat.test$p.value
  
  # stats
  test <- data %>% filter(gn %in% c("5x38", "7x24")) %>% 
    filter(sele == "cxcl12") %>%
    select(resno, gn, value) %>% unique()
  stat.test <- t.test(value ~ gn, data = test, exact = FALSE)
  stat.test$p.value
  
  
  ####
  data %>% filter(gn %in% c("5x38", "6x58",  "NTr.Cm1", "1x22")) %>%
    #filter(sub != "\\*") %>%
    select(resno, gn, sele, sub_mean) %>%
    filter(sele == "cxcl12") %>%
    unique() %>%
    ggplot(aes(gn, sub_mean)) +
    geom_boxplot() +
    theme_minimal()
  
  
  
  
