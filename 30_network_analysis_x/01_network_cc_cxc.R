# Name:     01_network_cc_cxc.R
# Updated:  20201103
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: ###############################
  
  # import, select those with evidence
  data <- read_csv("06_network/data/raw/ck_ckr_network.csv") %>%
    filter(evidence >=4) %>% select(ck, ckr) %>% unique() %>% 
    filter(ckr != "CXCR8") %>% filter(ckr != "H4 Receptor")
  
  # remove all but CC/CXC
  data <- data %>%
    filter(!grepl("ACKR", data$ckr)) %>% 
    filter(ckr != "XCR1") %>%
    filter(ckr != "CCRL2") %>%
    filter(ckr != "CX3CR1") %>%
    filter(ck != "XCL1") %>%
    filter(ck != "XCL2") %>%
    filter(ck != "CX3CL1")
    
  # add class labels
  data <- data %>% mutate(ck_class = case_when(
    grepl("CC", data$ck) ~ "cc",
    grepl("CXC", data$ck) ~ "cxc"
  )) %>% mutate(ckr_class = case_when(
    grepl("CC", data$ckr) ~ "cc",
    grepl("CXC", data$ckr) ~ "cxc"
  ))
  
  # count
  data <- data %>% select(ck_class, ckr_class) %>% count(ck_class, ckr_class)
  data <- data %>% mutate(fraction = n / (50+5+18))
  data <- data %>% unite(interaction, c(ck_class,ckr_class), sep = "_")
  
  # pie chart
  data %>%
    ggplot(aes(x = "", y= fraction, fill = interaction)) +
    geom_bar(width = 1,size = 1, stat="identity", color = "white") +
    coord_polar("y") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
  # stats
  test <- data.frame("CC" = c(50,0) , "CXC" = c(5,18))
  colnames(test) <- c("CC", "CXC")
  chisq.test(test)
  