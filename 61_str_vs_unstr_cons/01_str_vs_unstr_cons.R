# Name:     01_str_vs_unstr_cons.R
# Updated:  20210713
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

################################################################################

  data <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")

  # filter extreme N-term, batch approach (ie remove N-termini >10 residues)
  # & remove C-terminus
  data <- data %>% filter(!(ccn %in% c(
    "NTc.Cm11","NTc.Cm12","NTc.Cm13","NTc.Cm14","NTc.Cm15","NTc.Cm16","NTc.Cm17","NTc.Cm18","NTc.Cm19","NTc.Cm20","NTc.Cm21","NTc.Cm22","NTc.Cm23","NTc.Cm24","NTc.Cm25","NTc.Cm26","NTc.Cm27","NTc.Cm28","NTc.Cm29","NTc.Cm30","NTc.Cm31","NTc.Cm32","NTc.Cm33","NTc.Cm34","NTc.Cm35","NTc.Cm36","NTc.Cm37","NTc.Cm38","NTc.Cm39","NTc.Cm40","NTc.Cm41","NTc.Cm42","NTc.Cm43","NTc.Cm44","NTc.Cm45","NTc.Cm46","NTc.Cm47","NTc.Cm48","NTc.Cm49","NTc.Cm50","NTc.Cm51","NTc.Cm52","NTc.Cm53","NTc.Cm54","NTc.Cm55","NTc.Cm56","NTc.Cm57","NTc.Cm58","NTc.Cm59","NTc.Cm60","NTc.Cm61","NTc.Cm62","NTc.Cm63","NTc.Cm64","NTc.Cm65","NTc.Cm66","NTc.Cm67","NTc.Cm68","NTc.Cm69","NTc.Cm70","NTc.Cm71","NTc.Cm72","NTc.Cm73","NTc.Cm74","NTc.Cm75","NTc.Cm76","NTc.Cm77","NTc.Cm78","NTc.Cm79","NTc.Cm80","NTc.Cm81","NTc.Cm82","NTc.Cm83","NTc.Cm84","NTc.Cm85","NTc.Cm86","NTc.Cm87","NTc.Cm88","NTc.Cm89","NTc.Cm90","NTc.Cm91","NTc.Cm92","NTc.Cm93","NTc.Cm94","NTc.Cm95"
  ))) %>% filter(dom != "CT")  

  # select relevant cols
  # data <- data %>% select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # relabel N-term
  data <- data %>% mutate(nterm_or_not = case_when(
    dom == "NTc" ~ "nterm",
    dom != "NTc" ~ "core"
  ))
  
  # graph
  data %>%
    ggplot(aes(ortho_cons, all_cc_cxc_para)) +
    geom_jitter() +
    theme_minimal() +
    facet_grid(. ~ nterm_or_not)
  
  data %>%
    ggplot(aes(nterm_or_not, ortho_cons)) +
    geom_boxplot() +
    theme_minimal()
 
  data %>%
    filter(nterm_or_not == "nterm") %>%
    mean(ortho_cons)
  