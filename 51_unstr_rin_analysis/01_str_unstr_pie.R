# Name:     01_str_unstr_pie.R
# Updated:  20201125
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### STRUCTURED AND UNSTRUCTURED ANALYSIS - CHEMOKINE #########################

  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    # filter(file != "ngo" & file != "zheng") %>%
    select(source_gnccn, dom1, file)
    
  # add labels
  rin <- rin %>% mutate(str_ck = case_when(
    dom1 == "NTc" ~ "NTc",
    dom1 != "NTc" ~ "core_c"
  ))
  
  # summary stats
  rin.sum <- rin %>% count(file) # get no. interface positions per file
  colnames(rin.sum)[2] <- c("no_pos_per_file")
  str.unstr <- rin %>% count(str_ck, file)
  str.unstr <- left_join(str.unstr, rin.sum)
  rm(rin.sum, rin)
  str.unstr <- str.unstr %>% mutate(pct = n / no_pos_per_file)
  
  # violin/boxplot
  str.unstr %>%
    ggplot(aes(str_ck, pct)) +
    #geom_boxplot() +
    geom_violin(trim=FALSE) +
    #geom_point() +
    ylim(0,1) +
    theme_minimal() +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  
  # pie chart
  str.unstr <- str.unstr %>% group_by(str_ck) %>% summarise(mean = mean(pct), sd = sd(pct))
  
  str.unstr %>%
    ggplot(aes(x = "", y= mean)) +
    geom_bar(width = 1,size = 1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
##### STRUCTURED AND UNSTRUCTURED ANALYSIS - RECEPTOR ##########################
  
  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    # filter(file != "ngo" & file != "zheng") %>%
    select(target_gnccn, dom2, file)
  
  # add labels
  rin <- rin %>% mutate(str_ckr = case_when(
    dom2 == "NTr" ~ "NTr",
    dom2 == "ECL2" ~ "ECL2",
    dom2 != "NTr" & dom2 != "ECL2"  ~ "core_r"
  ))
  
  # summary stats
  rin.sum <- rin %>% count(file) # get no. interface positions per file
  colnames(rin.sum)[2] <- c("no_pos_per_file")
  str.unstr <- rin %>% count(str_ckr, file)
  str.unstr <- left_join(str.unstr, rin.sum)
  rm(rin.sum, rin)
  str.unstr <- str.unstr %>% mutate(pct = n / no_pos_per_file)
  
  # violin/boxplot
  str.unstr %>%
    ggplot(aes(str_ckr, pct)) +
    #geom_boxplot() +
    geom_violin(trim=FALSE) +
    #geom_point() +
    #ylim(0,1) +
    theme_minimal() +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2)
  
  # pie chart
  str.unstr <- str.unstr %>% group_by(str_ckr) %>% summarise(mean = mean(pct), sd = sd(pct))
  
  str.unstr %>%
    ggplot(aes(x = "", y= mean)) +
    geom_bar(width = 1,size = 1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  
  
  