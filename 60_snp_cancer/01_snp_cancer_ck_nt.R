# Name:     01_snp_cancer_nterm.R
# Updated:  20210524
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### 1: CHEMOKINE NTERM SNP ###################################################

  # import, reformat
  data <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")
  
  # filter extreme N-term, batch approach (ie remove N-termini >10 residues)
  # & remove C-terminus
  data <- data %>% filter(!(ccn %in% c(
    "NTc.Cm11","NTc.Cm12","NTc.Cm13","NTc.Cm14","NTc.Cm15","NTc.Cm16","NTc.Cm17","NTc.Cm18","NTc.Cm19","NTc.Cm20","NTc.Cm21","NTc.Cm22","NTc.Cm23","NTc.Cm24","NTc.Cm25","NTc.Cm26","NTc.Cm27","NTc.Cm28","NTc.Cm29","NTc.Cm30","NTc.Cm31","NTc.Cm32","NTc.Cm33","NTc.Cm34","NTc.Cm35","NTc.Cm36","NTc.Cm37","NTc.Cm38","NTc.Cm39","NTc.Cm40","NTc.Cm41","NTc.Cm42","NTc.Cm43","NTc.Cm44","NTc.Cm45","NTc.Cm46","NTc.Cm47","NTc.Cm48","NTc.Cm49","NTc.Cm50","NTc.Cm51","NTc.Cm52","NTc.Cm53","NTc.Cm54","NTc.Cm55","NTc.Cm56","NTc.Cm57","NTc.Cm58","NTc.Cm59","NTc.Cm60","NTc.Cm61","NTc.Cm62","NTc.Cm63","NTc.Cm64","NTc.Cm65","NTc.Cm66","NTc.Cm67","NTc.Cm68","NTc.Cm69","NTc.Cm70","NTc.Cm71","NTc.Cm72","NTc.Cm73","NTc.Cm74","NTc.Cm75","NTc.Cm76","NTc.Cm77","NTc.Cm78","NTc.Cm79","NTc.Cm80","NTc.Cm81","NTc.Cm82","NTc.Cm83","NTc.Cm84","NTc.Cm85","NTc.Cm86","NTc.Cm87","NTc.Cm88","NTc.Cm89","NTc.Cm90","NTc.Cm91","NTc.Cm92","NTc.Cm93","NTc.Cm94","NTc.Cm95"
  ))) %>% filter(dom != "CT")
  
  # select relevant cols
  data <- data %>% select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count)

  # filter high frequency alleles (30,000, ~10% frequency in population)
  data <- data %>% filter(snp_freq_count < 30000)
  
  # relabel N-term
  data <- data %>% mutate(nterm_or_not = case_when(
    dom == "NTc" ~ "nterm",
    dom != "NTc" ~ "core"
  ))

  # setup chi squared contingency table
  test <- data %>% group_by(nterm_or_not) %>% 
    filter(!is.na(snp_freq_count)) %>% 
    summarise(n=dplyr::n(), snp_freq_count = sum(snp_freq_count)) %>% 
    ungroup()
  chisq.test(c(20737, 66104), p=c(427, 2488)/(427 + 2488))

  # define function to graph expected and observed chi square values for both groups
  ChiToGraph <- function(DF, TIER, NON){
    data <- DF %>% gather(key, value, 2:3)
    data <- data %>% dplyr::group_by(key) %>% dplyr::mutate(total = sum(value)) %>% dplyr::ungroup()
    data <- data %>% mutate(fraction = value / total)
    data <- data.frame("EO" = c("O", "E", "O", "E"), 
                       "tier" = c(TIER, TIER, NON, NON),
                       "value" = c(data$value[4], data$fraction[2]*data$total[3], data$value[3], data$fraction[1]*data$total[3]) )
    return(data)
    rm(data)
  }
  
  # PLOT - EXPECTED VS. OBSERVED
  test <- ChiToGraph(test, "nterm", "core")
  
  order = c("nterm", "core")
  order2 = c("E", "O")
  test$tier <- factor(test$tier, levels = order)
  test$EO <- factor(test$EO, levels = order2)
  test %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # Hypothesis: chemokines have more SNPs than expected in N-terminus
  # Results: True! When you remove C-terminal region

  # (best would be to then show unique activity for Nterm SNP/cancer)
  # coulds show variation in Nterm (ideally motif, ie conserved) leads to
  # functional innovation in human DISEASE
  
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(-0.3,1)
  
##### 2: CHEMOKINE NTERM CANCER ################################################
  
  # import, reformat
  data <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")
  
  # filter extreme N-term, batch approach (ie remove N-termini >10 residues)
  # & remove C-terminus
  data <- data %>% filter(!(ccn %in% c(
    "NTc.Cm11","NTc.Cm12","NTc.Cm13","NTc.Cm14","NTc.Cm15","NTc.Cm16","NTc.Cm17","NTc.Cm18","NTc.Cm19","NTc.Cm20","NTc.Cm21","NTc.Cm22","NTc.Cm23","NTc.Cm24","NTc.Cm25","NTc.Cm26","NTc.Cm27","NTc.Cm28","NTc.Cm29","NTc.Cm30","NTc.Cm31","NTc.Cm32","NTc.Cm33","NTc.Cm34","NTc.Cm35","NTc.Cm36","NTc.Cm37","NTc.Cm38","NTc.Cm39","NTc.Cm40","NTc.Cm41","NTc.Cm42","NTc.Cm43","NTc.Cm44","NTc.Cm45","NTc.Cm46","NTc.Cm47","NTc.Cm48","NTc.Cm49","NTc.Cm50","NTc.Cm51","NTc.Cm52","NTc.Cm53","NTc.Cm54","NTc.Cm55","NTc.Cm56","NTc.Cm57","NTc.Cm58","NTc.Cm59","NTc.Cm60","NTc.Cm61","NTc.Cm62","NTc.Cm63","NTc.Cm64","NTc.Cm65","NTc.Cm66","NTc.Cm67","NTc.Cm68","NTc.Cm69","NTc.Cm70","NTc.Cm71","NTc.Cm72","NTc.Cm73","NTc.Cm74","NTc.Cm75","NTc.Cm76","NTc.Cm77","NTc.Cm78","NTc.Cm79","NTc.Cm80","NTc.Cm81","NTc.Cm82","NTc.Cm83","NTc.Cm84","NTc.Cm85","NTc.Cm86","NTc.Cm87","NTc.Cm88","NTc.Cm89","NTc.Cm90","NTc.Cm91","NTc.Cm92","NTc.Cm93","NTc.Cm94","NTc.Cm95"
  ))) %>% filter(dom != "CT")
  
  # select relevant cols
  data <- data %>% select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # relabel N-term
  data <- data %>% mutate(nterm_or_not = case_when(
    dom == "NTc" ~ "nterm",
    dom != "NTc" ~ "core"
  ))
  
  # setup chi squared contingency table
  test <- data %>% group_by(nterm_or_not) %>% 
    filter(!is.na(cancer_mut_count)) %>% 
    summarise(n=dplyr::n(), cancer_mut_count = sum(cancer_mut_count)) %>% 
    ungroup()
  chisq.test(c(65, 446), p=c(428, 2491)/(428 + 2491))

    
  # define function to graph expected and observed chi square values for both groups
  ChiToGraph <- function(DF, TIER, NON){
    data <- DF %>% gather(key, value, 2:3)
    data <- data %>% dplyr::group_by(key) %>% dplyr::mutate(total = sum(value)) %>% dplyr::ungroup()
    data <- data %>% mutate(fraction = value / total)
    data <- data.frame("EO" = c("O", "E", "O", "E"), 
                       "tier" = c(TIER, TIER, NON, NON),
                       "value" = c(data$value[4], data$fraction[2]*data$total[3], data$value[3], data$fraction[1]*data$total[3]) )
    return(data)
    rm(data)
  }
  
  # graph results
  test <- ChiToGraph(test, "nterm", "core")
  
  order = c("nterm", "core")
  order2 = c("E", "O")
  test$tier <- factor(test$tier, levels = order)
  test$EO <- factor(test$EO, levels = order2)
  test %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(-0.3,1)
  
  
  