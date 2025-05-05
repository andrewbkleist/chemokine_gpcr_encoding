# Name:     04_core_vs_nterm_seq_sim.R
# Updated:  20201214
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggrepel)

##### 1: CORE VS. NTERM SEQ SIM - CHEMOKINE  ###################################
  
  # import
  data <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")
  
  # add labels
  data <- data %>% 
    mutate(str_ck = case_when(
      dom == "NTc" ~ "NTc",
      dom != "NTc" ~ "core_c"
    )) 
  
  # subset
  nterm <- data %>% filter(str_ck == "NTc")
  core <- data %>% filter(str_ck != "NTc") %>% filter(dom != "CT")
  data <- bind_rows(nterm, core)
  rm(nterm, core)
  
  data <- data %>% filter(ortho_cons != 0)
  
  # plot
  data %>%
    ggplot(aes(str_ck, ortho_cons)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot() +
    #geom_jitter(shape = 21, colour = "black", fill = "white", size = 2, stroke = 0.5) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, colour = "black", fill = "white", stroke = 0.5) +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    #geom_boxplot() +
    ylim(0,1) +
    theme_minimal() 

  # stats  
  sum <- wilcox.test(ortho_cons ~ str_ck, data = data, exact = FALSE)
  sum$p.value
  
  
  
##### 2: CORE VS. NTERM SEQ SIM - RECEPTOR  ####################################
  
  # import
  data <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv")
  
  # add labels
  data <- data %>% 
    mutate(str_ckr = case_when(
      dom == "NTr" ~ "NTr",
      dom == "ECL2" ~ "ECL2",
      dom != "NTr" & dom != "ECL2"  ~ "core_r"
    ))
  
  # subset
  nterm <- data %>% filter(str_ckr == "NTc")
  core <- data %>% filter(str_ckr != "NTc") %>% filter(dom != "CT")
  data <- bind_rows(nterm, core)
  rm(nterm, core)
  
  data <- data %>% filter(ortho_cons != 0)
  
  # plot
  data %>%
    ggplot(aes(str_ckr, ortho_cons)) +
    #geom_violin(trim = FALSE) +
    geom_boxplot() +
    #geom_jitter(shape = 21, colour = "black", fill = "white", size = 2, stroke = 0.5) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, colour = "black", fill = "white", stroke = 0.5) +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    #geom_boxplot() +
    ylim(0,1) +
    theme_minimal() 
  
  # stats  
  
  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  core.nt <- data %>% filter(str_ckr != "ECL2")
  core.nt <- wilcox.test(ortho_cons ~ str_ckr, data = core.nt, exact = FALSE)
  core.nt$p.value
  
  core.ecl2 <-  data %>% filter(str_ckr != "NTr")
  core.ecl2 <- wilcox.test(ortho_cons ~ str_ckr, data = core.ecl2, exact = FALSE)
  core.ecl2$p.value
  
  nt.ecl2 <-  data %>% filter(str_ckr != "core_r")
  nt.ecl2 <- wilcox.test(ortho_cons ~ str_ckr, data = nt.ecl2, exact = FALSE)
  nt.ecl2$p.value
  

  