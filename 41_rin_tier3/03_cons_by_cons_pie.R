# Name:     03_cons_by_cons_pie.R
# Updated:  20210415
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: CHEMOKINE PIE CHARTS ##################################################
  
  # import
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # select and pie chart - CHEMOKINE T1
  rin %>% 
    filter(file == "5uiw") %>% 
    filter(source_gnccn %in% c("b1b2.9")) %>%
    select(source_gnccn, all_para_ck) %>% 
    unique() %>%
    ggplot(aes(x = "", y= all_para_ck)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ source_gnccn)
  
  # select and pie chart - CHEMOKINE T2
  # note that you must scale T2 since it is on a 0.5-1 scale
  rin %>% 
    filter(file == "5uiw") %>% 
    filter(source_gnccn %in% c("b1b2.9")) %>%
    select(source_gnccn, cc_cxc_lr_ck) %>% 
    unique() %>%
    mutate(cc_cxc_lr_ck = cc_cxc_lr_ck - 0.5) %>%
    ggplot(aes(x = "", y= cc_cxc_lr_ck)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0, .5) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ source_gnccn)
  
  # select and pie chart - CHEMOKINE T3
  rin %>% 
    filter(file == "5uiw") %>% 
    filter(source_gnccn %in% c("b1b2.9")) %>%
    select(source_gnccn, ortho_cons_ck) %>% 
    unique() %>%
    ggplot(aes(x = "", y= ortho_cons_ck)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ source_gnccn)
  
  # select and pie chart - CHEMOKINE T3
  rin %>% 
    filter(file == "6wwz") %>% 
    filter(source_gnccn %in% c("b1b2.9")) %>%
    select(source_gnccn, ortho_cons_ck) %>% 
    unique() %>%
    ggplot(aes(x = "", y= ortho_cons_ck)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ source_gnccn)
  
  
  
##### 2: RECEPTOR PIE CHARTS ###################################################
  
  # import
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # select and pie chart -  T1
  rin %>% 
    filter(file == "5uiw") %>% 
    filter(target_gnccn %in% c("ECL2.Cp4")) %>%
    select(target_gnccn, all_para_ckr) %>% 
    unique() %>%
    ggplot(aes(x = "", y= all_para_ckr)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ target_gnccn)
  
  # select and pie chart -  T1
  rin %>% 
    filter(file == "6wwz") %>% 
    filter(target_gnccn %in% c("ECL2.Cp3")) %>%
    select(target_gnccn, all_para_ckr) %>% 
    unique() %>%
    ggplot(aes(x = "", y= all_para_ckr)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ target_gnccn)
  
  # select and pie chart -  T2
  # note that you must scale T2 since it is on a 0.5-1 scale
  # rin %>% 
  #   filter(file == "5uiw") %>% 
  #   filter(target_gnccn %in% c("ECL2.Cp4")) %>%
  #   select(target_gnccn, cc_cxc_lr_ckr) %>% 
  #   unique() %>%
  #   mutate(cc_cxc_lr_ckr = cc_cxc_lr_ckr - 0.5) %>%
  #   ggplot(aes(x = "", y= cc_cxc_lr_ckr)) +
  #   geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  #   coord_polar("y") +
  #   ylim(0, .5) +
  #   theme_classic() +
  #   theme(axis.line = element_blank(),
  #         axis.text = element_blank(),
  #         axis.ticks = element_blank()) +
  #   facet_grid(. ~ target_gnccn)
  
  # select and pie chart -  T3
  rin %>% 
    filter(file == "5uiw") %>% 
    filter(target_gnccn %in% c("ECL2.Cp4")) %>%
    select(target_gnccn, ortho_cons_ckr) %>% 
    unique() %>%
    ggplot(aes(x = "", y= ortho_cons_ckr)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ target_gnccn)
  
  
  # select and pie chart -  T3
  rin %>% 
    filter(file == "6wwz") %>% 
    filter(target_gnccn %in% c("ECL2.Cp3")) %>%
    select(target_gnccn, ortho_cons_ckr) %>% 
    unique() %>%
    ggplot(aes(x = "", y= ortho_cons_ckr)) +
    geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
    coord_polar("y") +
    ylim(0,1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank()) +
    facet_grid(. ~ target_gnccn)
  
  
  
  