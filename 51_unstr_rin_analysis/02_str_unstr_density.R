# Name:     02_str_unstr_density.R
# Updated:  20201125
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### STRUCTURED AND UNSTRUCTURED ANALYSIS - CHEMOKINE #########################
  
  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(file != "6meo") %>%
    select(source_gnccn, dom1, file)
  
  # add labels
  rin <- rin %>% mutate(str_ck = case_when(
    dom1 == "NTc" ~ "NTc",
    dom1 != "NTc" ~ "core_c"
  ))
  
  # summary stats
  contacts <- rin %>% count(source_gnccn, file, str_ck)

  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  sum <- wilcox.test(n ~ str_ck, data = contacts, exact = FALSE)
  sum$p.value
  
  # plot
  contacts %>%
    ggplot(aes(str_ck, n)) +
    geom_boxplot() +
    #geom_violin(trim = FALSE) +
    ylim(0,7) +
    theme_minimal() 
    #stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  
  
##### STRUCTURED AND UNSTRUCTURED ANALYSIS - RECEPTOR #########################
  
  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(file != "6meo") %>%
    select(target_gnccn, dom2, file)
  
  # # add labels
  # rin <- rin %>% mutate(str_ckr = case_when(
  #   dom2 == "NTr" ~ "NTr",
  #   dom2 == "ECL2" ~ "ECL2",
  #   dom2 != "NTr" & dom2 != "ECL2"  ~ "core_r"
  # ))
  
  # add labels
  rin <- rin %>% mutate(str_ckr = case_when(
    dom2 == "NTr" ~ "NTr",
    dom2 != "NTr" ~ "core_r"
  ))
  
  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  contacts <- rin %>% count(target_gnccn, file, str_ckr)
  
  ckr.nt_core <- filter(contacts, str_ckr == "NTr" | str_ckr == "core_r" )
  ckr.nt_core <- wilcox.test(n ~ str_ckr, data = ckr.nt_core, exact = FALSE)
  ckr.nt_core$p.value # significant
  
  # ckr.nt_ecl2 <- filter(contacts, str_ckr == "NTr" | str_ckr == "ECL2" )
  # ckr.nt_ecl2 <- wilcox.test(n ~ str_ckr, data = ckr.nt_ecl2, exact = FALSE)
  # ckr.nt_ecl2$p.value # not significant
  # 
  # ckr.core_ecl2 <- filter(contacts, str_ckr == "core_r" | str_ckr == "ECL2" )
  # ckr.core_ecl2 <- wilcox.test(n ~ str_ckr, data = ckr.core_ecl2, exact = FALSE)
  # ckr.core_ecl2$p.value # significant
  
  # plot
  contacts %>%
    ggplot(aes(str_ckr, n)) +
    geom_boxplot() +
    #geom_violin(trim = FALSE) +
    ylim(0,7) +
    theme_minimal() 
    #stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  
  
