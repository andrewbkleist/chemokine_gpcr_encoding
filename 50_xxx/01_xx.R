# Name:     get_seq_struct_differences.R
# Updated:  20230206
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##############################################################################





data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(file %in% c("5uiw","7o7f", "7f1r", "zheng", "7vl9", "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "7sk3", "7xbx")) %>%
  # select all complexes containing only human CK/CKR (12 complexes)
  select(file, source_gnccn, target_gnccn, dom1, dom2, all_para_ck, all_para_ckr) %>%
  unique()
data$mean <- ave(data$all_para_ck, data$all_para_ckr)
get.pdb <- data %>% count(source_gnccn, target_gnccn)
data <- left_join(data, get.pdb)

data %>%
  ggplot(aes(n, mean)) +
  geom_jitter(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5, width = 0.4, height = 0.05) +
  scale_x_continuous(limits=c(0, 9), breaks=c(2,4,6,8,10)) +
  scale_y_continuous(limits=c(0, 1)) +
  theme_minimal()

data <- data %>%
  mutate(str_unstr_ck = case_when(
    dom1 == "NTc" ~ "NTc",
    dom1 != "NTc" ~ "str"
  )) %>%
  mutate(str_unstr_ckr = case_when(
    dom2 == "NTr" ~ "NTr",
    dom2 == "ECL2" ~ "ECL2",
    dom2 != c("NTc", "ECL2")  ~ "str"
  ))

# filter and count
data <- data %>% filter(n == 1) %>% filter(mean < 0.5) %>% count(str_unstr_ck, str_unstr_ckr)
data <- data %>% mutate(pct = n / (sum(n)))

data$n <- factor(data$n, levels = data$n[order(data$n)])
data %>%
ggplot(aes(x = "", y= pct)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

#----

GetContactVenn <- function(PDB1, PDB2){
  # Defines percentage of contacts shared (intersection) for 2 PDBs
  # compared to unique contacts in each PDB
  set1 <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file == PDB1) %>%
    select(file, source_gnccn, target_gnccn) %>%
    unique() %>%
    select(-file)
  
  set2 <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
    filter(file == PDB2) %>%
    select(file, source_gnccn, target_gnccn) %>%
    unique() %>%
    select(-file)
  
  set1.not.set2 <- setdiff(set1, set2)
  set1.not.set2$type <- c(paste0(PDB1, "_not_", PDB2))
  set2.not.set1 <- setdiff(set2, set1)
  set2.not.set1$type <- c(paste0(PDB2, "_not_", PDB1))
  set1.i.set2 <- intersect(set1, set2)
  set1.i.set2$type <- c(paste0(PDB1, "_", PDB2))
  master <- rbind(set1.not.set2, set2.not.set1, set1.i.set2)
  rm(set1, set2)
  
  return(master)
}

test <- GetContactVenn("5uiw", "7vl9")
test <- GetContactVenn("5uiw", "7o7f")
test <- GetContactVenn("5wb2", "4xt1")



#--

data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(file %in% c("5uiw","7vl9", "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "7sk3", "7xbx")) %>%
  select(file, no_pdb, source_gnccn, target_gnccn, all_para_ck, all_para_ckr) %>%
  unique()

data$mean <- ave(data$all_para_ck, data$all_para_ckr)


data %>%
  ggplot(aes(no_pdb, mean)) +
  geom_jitter(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5, width = 0.5, height = 0.05) +
  ylim(0,1) +
  xlim(1,10) +
  theme_minimal()
