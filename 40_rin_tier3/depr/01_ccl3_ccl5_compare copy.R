# Name:     01_ccl3_ccl5_compare.R
# Updated:  20230123
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### FUNCTIONS ################################################################

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
  set2.not.set1 <- setdiff(set2, set1)
  set1.i.set2 <- intersect(set1, set2)
  rm(set1, set2)
  
  return(nrow(set1.i.set2)/(nrow(set1.not.set2) + nrow(set2.not.set1) + nrow(set1.i.set2)))
}

################################################################################

# define all PDBs
pdbs <- c("5uiw", "7o7f", "7f1r", "zheng",
          "7vl9", "7xa3", "7f1t", "6wwz",
          "6lfo", "ngo",
          "7sk3", "7xbx",
          "4rws", "4xt1", "5wb2")

# create empty df and fill with pairwise RIN PIDs
df <- NULL

for (i in pdbs){
  for (j in pdbs){
    data <- GetContactVenn(i, j)
    temp <- data.frame(data, paste0(i), paste0(j))
    df <- rbind(df, temp)
  }
}
rm(temp, i, j, data)
colnames(df) <- c("rin_pid", "pdb1", "pdb2")

# add chemokine, receptor information
chemokine <- c("CCL5", "CCL5", "CCL5", "CCL5", 
               "CCL15", "CCL2", "CCL3", "CCL20", 
               "CXCL8", "CXCL12", "CXCL12",
               "CX3CL1", "vMIPII", "CX3CL1", "CX3CL1"
               )

receptor <- c("CCR5", "CCR5", "CCR5", "CCR5", 
              "CCR1", "CCR2", "CCR5", "CCR6",
              "CXCR2", "CXCR4",
              "ACKR3", "CX3CR1",
              "CXCR4", "US28", "US28")
lookup <- data.frame(pdbs, chemokine, receptor)

# map chemokine/receptor names to df, unite columns
df$ck1 <- lookup$chemokine[match(unlist(df$pdb1), lookup$pdbs)]
df$ckr1 <- lookup$receptor[match(unlist(df$pdb1), lookup$pdbs)]
df$ck2 <- lookup$chemokine[match(unlist(df$pdb2), lookup$pdbs)]
df$ckr2 <- lookup$receptor[match(unlist(df$pdb2), lookup$pdbs)]
df <- df %>% unite(ck1_ck2, c(ck1, ck2))
df <- df %>% unite(ckr1_ckr2, c(ckr1, ckr2))

# import distances, unite columns
ck.dist <- read_csv("06_network/output/network_similarity_scores_ck.csv")
ckr.dist <- read_csv("06_network/output/network_similarity_scores_ckr.csv")

ck.dist <- ck.dist %>% unite(ck1_ck2, c(a, b))
ckr.dist <- ckr.dist %>% unite(ckr1_ckr2, c(a, b))

# map distances to df
df$ck_sim <- ck.dist$sim_nl[match(unlist(df$ck1_ck2), ck.dist$ck1_ck2)]
df$ckr_sim <- ckr.dist$sim_nl[match(unlist(df$ckr1_ckr2), ckr.dist$ckr1_ckr2)]

# clean up df
df <- df %>% filter(pdb1 !=pdb2)
df <- df %>% filter(!(pdb1 %in% c("4rws", "zheng", "7o7f", "7f1r", "4xt1", "5wb2")))
df <- df %>% filter(!(pdb2 %in% c("4rws", "zheng", "7o7f", "7f1r", "4xt1", "5wb2")))


# test
df$arb <- c(1)
df <- df %>% select(rin_pid:pdb2)
df %>%
  as.matrix %>%
  cor %>% {(function(x){x[upper.tri(x)]<-NA; x})(.)} %>%
  as.data.frame %>%
  rownames_to_column(var = 'var1') %>%
  gather(var2, value, -var1) %>%
  filter(value > 0.8 | value < -0.8) %>%
  filter(value != 1)

x <- cor(df)
x[ upper.tri(x, diag = TRUE) | abs(x) < 0.8  ] <- NA
na.omit(data.frame(as.table(x)))

df %>%
  ggplot(aes(rin_pid, arb)) +
  geom_jitter(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  xlim(0, 0.5) +
  geom_vline(xintercept = c(0.1, 0.2)) +
  theme_minimal()

df %>%
  ggplot(aes(x=rin_pid)) +
geom_histogram(binwidth = 0.01) +
  xlim(0, 0.5) +
  theme_minimal()


# plot
df %>%
  # filter(pdb1 != pdb2) %>%
  ggplot(aes(ckr_sim, ck_sim, size = rin_pid)) +
  geom_point(shape = 21, colour = "black", fill = "white", stroke = 0.5) +
  scale_size(range = c(1, 20)) +
  xlim(0,1) +
  ylim(0,1) +
  # geom_text() +
  theme_minimal()

df %>%
  # filter(pdb1 != pdb2) %>%
  ggplot(aes(rin_pid, ckr_sim)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  # scale_size(range = c(1, 15)) +
  # xlim(0,1) +
  # ylim(0,1) +
  theme_minimal()

df %>%
  # filter(pdb1 != pdb2) %>%
  ggplot(aes(rin_pid, ck_sim)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  # scale_size(range = c(1, 15)) +
  # xlim(0,1) +
  # ylim(0,1) +
  theme_minimal()




+++++
a <- GetContactVenn("7o7f", "7vl9")
b <- GetContactVenn("7o7f", "7xa3")
c <- GetContactVenn("7o7f", "7f1t")
d <- GetContactVenn("7o7f", "6wwz")
e <- GetContactVenn("7vl9", "7xa3")
f <- GetContactVenn("7vl9", "7f1t")
g <- GetContactVenn("7vl9", "6wwz")
h <- GetContactVenn("7xa3", "7f1t")
i <- GetContactVenn("7xa3", "6wwz")
j <- GetContactVenn("7f1t", "6wwz")

test <- GetContactVenn("6wwz", "6lfo")

test <- GetContactVenn("4xt1", "5wb2")












  
  
rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
  # filter(file %in% c("7f1t", "5uiw", "6wwz")) %>% 
  select(file, source_gnccn, target_gnccn) %>%
  unique()

set1 <- rin %>% filter(file == "7f1t") %>% select(-file) # CCL3:CCR5
set2 <- rin %>% filter(file == "7vl9") %>% select(-file) # CCL5[5p7]: CCR5

set1.not.set2 <- setdiff(set1, set2)
set2.not.set1 <- setdiff(set2, set1)
set1.i.set2 <- intersect(set1, set2)






  
  

