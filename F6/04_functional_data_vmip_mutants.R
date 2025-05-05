# Name:     06_functional_data_vmip_mutants.R
# Updated:  20230109
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################
  
  # import data
  data <- read_csv("35_viral_ackr/data/processed/vmipii_data_summary.csv")
  
  # subset data by receptor, compare to WT vMIP-II state
  ccr3 <- data %>% filter(receptor=="ccr3") %>% filter(metric != "Emax")
  ccr3 <- ccr3 %>% mutate(fold_change = log(1.35e-09/value))
  
  ccr5 <- data %>% filter(receptor=="ccr5") %>% filter(metric != "Emin")
  ccr5 <- ccr5 %>% mutate(fold_change = log(1.165e-08/value))
  
  cxcr4 <- data %>% filter(receptor=="cxcr4") %>% filter(metric != "Emin")
  cxcr4 <- cxcr4 %>% mutate(fold_change = log(9.34e-09/value))
  
  # recombine
  data <- rbind(ccr3, ccr5, cxcr4)
  rm(ccr3, ccr5, cxcr4)
  
  # plot
  order <- c("WT", "R7I", "K10T", "L13F", "R7I_K10T_L13F")
  data$mutant <- factor(data$mutant, levels = order)
  
  
  # by mutant circular
  # data %>%
  #   filter(mutant != "WT") %>%
  #   ggplot(aes(mutant, fold_change, color = receptor)) +
  #   geom_point() +
  #   ylim(-3, 2) +
  #   # geom_line(group = 1) +
  #   coord_polar() +
  #   facet_grid(. ~ receptor) +
  #   theme_minimal()
  
  # by mutant plain
  data %>%
    filter(mutant != "WT") %>%
    ggplot(aes(mutant, fold_change, color = receptor)) +
    geom_point(shape = 21, fill = "white", size = 5, stroke = 0.5) +
    # ylim(-3, 1) +
    # geom_line(group = 1) +
    # coord_polar() +
    # facet_grid(. ~ receptor) +
    theme_minimal()
  
  # data %>%
  #   # filter(mutant != "WT") %>%
  #   ggplot(aes(receptor, fold_change, color = mutant)) +
  #   geom_point(shape = 21, fill = "white", size = 5, stroke = 0.5) +
  #   ylim(-3, 1) +
  #   # geom_line(group = 1) +
  #   coord_polar() +
  #   # facet_grid(. ~ receptor) +
  #   theme_minimal()
  
  data %>%
    filter(mutant != "WT") %>%
    ggplot(aes(mutant, fold_change)) +
    # geom_bar(stat = "identity") +
    geom_segment(aes(x=mutant, xend=mutant, y=0, yend=fold_change)) +
    geom_point(shape = 21, fill = "white", size = 5, stroke = 0.5) +
    # ylim(-3, 1) +
    # geom_segment(aes(x = receptor, y = 0, xend = receptor, yend = fold_change)) +
    # coord_polar() +
    facet_grid(. ~ receptor) +
    theme_minimal()
  
  ggsave(filename = "mutant_fold_ec50.pdf", 
         plot = last_plot(), path = "F6/output/",
         width = 5,
         height = 3)



################################################################################

  # import data
  data <- read_csv("35_viral_ackr/data/processed/vmipii_data_summary.csv")
  
  # subset data by receptor, compare to WT vMIP-II state
  ccr3 <- data %>% filter(receptor=="ccr3") %>% filter(metric != "Emax")
  ccr3 <- ccr3 %>% mutate(fold_change = log(1.35e-09/value))
  
  ccr5 <- data %>% filter(receptor=="ccr5") %>% filter(metric != "Emin")
  ccr5 <- ccr5 %>% mutate(fold_change = log(1.165e-08/value))
  
  cxcr4 <- data %>% filter(receptor=="cxcr4") %>% filter(metric != "Emin")
  cxcr4 <- cxcr4 %>% mutate(fold_change = log(9.34e-09/value))
  
  # recombine
  data1 <- rbind(ccr3, ccr5, cxcr4)
  data1$type <- c("ec50")
  rm(ccr3, ccr5, cxcr4)
  
  # subset data by receptor, compare to WT vMIP-II state
  ccr3 <- data %>% filter(receptor=="ccr3") %>% filter(metric != "EC50")
  ccr3 <- ccr3 %>% mutate(fold_change = log(abs(value/1.206)))
  
  ccr5 <- data %>% filter(receptor=="ccr5") %>% filter(metric != "EC50")
  ccr5 <- ccr5 %>% mutate(fold_change = log(abs(value/-0.1273)))
  
  cxcr4 <- data %>% filter(receptor=="cxcr4") %>% filter(metric != "EC50")
  cxcr4 <- cxcr4 %>% mutate(fold_change = log(abs(value/-18.39)))
  
  # recombine
  data2 <- rbind(ccr3, ccr5, cxcr4)
  data2$type <- c("eminmax")
  rm(ccr3, ccr5, cxcr4)
  
  # combine
  data <- rbind(data1, data2)
  data <- data  %>% select(-metric, -value)
  data <- data %>% pivot_wider(names_from = type, values_from = fold_change)
  
  # plot
  data %>%
    ggplot(aes(ec50, eminmax, color = receptor)) +
    geom_point(shape = 21, fill = "white", size = 5, stroke = 0.5) +
    theme_minimal()
    
  
  
  
  

  
  
  