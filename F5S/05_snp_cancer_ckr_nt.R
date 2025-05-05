# Name:     02_snp_cancer_ckr_nt.R
# Updated:  20210705
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### 1: RECEPTOR NTERM SNP ####################################################
  
  # import, reformat
  data <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv")
  
  # filter out extreme C-terminus
  data <- data %>% filter(dom != "CT")
  
  # filter ECL2 (considered  separately)
  # data <- data %>% filter(dom != "ECL2")
  
  # select relevant cols
  data <- data %>% select(protein, gn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # filter high frequency alleles (30,000, ~10% frequency in population)
  data <- data %>% filter(snp_freq_count < 30000)
  
  # relabel N-term
  data <- data %>% mutate(nterm_or_not = case_when(
    dom == "NTr" ~ "nterm",
    dom != "NTr" ~ "core"
  ))
  
  # setup chi squared contingency table
  test <- data %>% group_by(nterm_or_not) %>% 
    filter(!is.na(snp_freq_count)) %>% 
    summarise(n=dplyr::n(), snp_freq_count = sum(snp_freq_count)) %>% 
    ungroup()
  chisq.test(c(27078, 134932), p=c(708, 6737)/(708 + 6737))
  
  # count by dom - for reference only
  # no_per_dom <- data %>% group_by(dom) %>% summarize(count = count(dom), no_per_dom = sum(snp_freq_count)) %>%
  #   ungroup()
  
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
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(-0.3,1)
  
  ggsave(filename = "ckr_snp.pdf", 
         plot = last_plot(), path = "F5S/output/",
         width = 3,
         height = 6)
  
  
##### 2: RECEPTOR NTERM CANCER #################################################
  
  # import, reformat
  data <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv")
  
  # filter out extreme C-terminus
  data <- data %>% filter(dom != "CT")
  
  # filter ECL2 (considered  separately)
  data <- data %>% filter(dom != "ECL2")
  
  # select relevant cols
  data <- data %>% select(protein, gn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # relabel N-term
  data <- data %>% mutate(nterm_or_not = case_when(
    dom == "NTr" ~ "nterm",
    dom != "NTr" ~ "core"
  ))
  
  # setup chi squared contingency table
  test <- data %>% group_by(nterm_or_not) %>% 
    filter(!is.na(cancer_mut_count)) %>% 
    summarise(n=dplyr::n(), cancer_mut_count = sum(cancer_mut_count)) %>% 
    ungroup()
  chisq.test(c(140, 1095), p=c(709, 6323)/(709 + 6323))
  
  # count by dom - for reference only
  # no_per_dom <- data %>% group_by(dom) %>% summarize(count = count(dom), no_per_dom = sum(snp_freq_count)) %>%
  #   ungroup()
  
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
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(-0.3,1)
  
  
  
  