# Name:     01_ck_interface_t2_score_ridge.R
# Updated:  20230106
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggridges)


##### 1: SCORE COMPARISONS INTERFACE - RIDGE PLOT ##############################
  
  # plot spectrum of T2-T2 scores from complexes
  pdbs <- c("7vl9", "7xa3", "7f1t", "6wwz", "zheng", "6lfo", "ngo", "4rws")
    
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(cc_cxc_lr_ck >= 0.75) %>%
    filter(file %in% pdbs) %>%
    select(ck, source_gnccn, cc_cxc_lr_score_ck) %>% 
    unique()
  
  order <- c("ccl2", "ccl3", "ccl5", "ccl15", "ccl20","vmipii",  "cxcl12", "cxcl8")
  rin$ck <- factor(rin$ck, levels = rev(order))
  
  rin %>%
    ggplot(aes(x = cc_cxc_lr_score_ck, y = ck)) +
    geom_density_ridges(fill = "steelblue3") +
    theme_minimal() 
  

##### 2: SCORE COMPARISONS ALL POSITIONS - RIDGE PLOT ##########################
  
data <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
vmipii <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
data <- rbind(data, vmipii)
t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  filter(mean >= 0.75)
t2 <- as.character(t2$motif)
rm(vmipii)
data %>%
  filter(position %in% t2) %>%
  filter(protein %in% c("ccl2", "ccl3", "ccl5", "ccl15", "ccl20","vmip2xhhv8p",  "cxcl12", "cxcl8")) %>%
  # filter(position %in% c("NTc.Cm10","NTc.Cm9","NTc.Cm8","NTc.Cm7","NTc.Cm6","NTc.Cm5","NTc.Cm4","NTc.Cm3","NTc.Cm2","NTc.Cm1", "NTc.Cm0", "CX.1", "CX.2", "CX.3", "CX.4", "CX.5","cxb1.1","cxb1.2","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.10","cxb1.11","cxb1.14","cxb1.15","cxb1.16","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.4","b1b2.6","b1b2.7","b1b2.8", "b1b2.9","b1b2.10","b1b2.12","b1b2.13","b1b2.14","b1b2.16","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.2","b2b3.3","b2b3.4","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10")) %>%
  ggplot(aes(x = score, y = protein)) +
  geom_density_ridges(fill = "steelblue3") +
  theme_minimal() 

test %>%
  filter(protein %in% c("vmip2xhhv8p")) %>%
  filter(position %in% c("NTc.Cm10","NTc.Cm9","NTc.Cm8","NTc.Cm7","NTc.Cm6","NTc.Cm5","NTc.Cm4","NTc.Cm3","NTc.Cm2","NTc.Cm1", "NTc.Cm0", "CX.1", "CX.2", "CX.3", "CX.4", "CX.5","cxb1.1","cxb1.2","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.10","cxb1.11","cxb1.14","cxb1.15","cxb1.16","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.4","b1b2.6","b1b2.7","b1b2.8", "b1b2.9","b1b2.10","b1b2.12","b1b2.13","b1b2.14","b1b2.16","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.2","b2b3.3","b2b3.4","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10")) %>%
  ggplot(aes(x = score, y = protein)) +
  geom_density_ridges(fill = "steelblue3") +
  theme_minimal() 
