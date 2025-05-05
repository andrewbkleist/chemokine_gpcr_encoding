# Updated:  20230123
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
#   filter(file != "6meo") %>%
#   select(source_gnccn, target_gnccn, dom1, dom2) %>%
#   unique()

data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(no_pdb == 1) %>%
  filter(all_para_ck < 0.5) %>%
  filter(all_para_ckr < 0.5) %>%
  filter(file != "6meo") %>%
  filter(file != "ngo") %>%
  
  select(source_gnccn, target_gnccn, dom1, dom2) %>%
  unique() # 176 contacts (of 478 unique contacts = 37%)

data <- data %>%
  mutate(str_unstr_ck = case_when(
    dom1 == "NTc" ~ "unstr",
    dom1 != "NTc" ~ "str",
  )) %>% mutate(str_unstr_ckr = case_when(
    dom2 == "NTr" ~ "unstr",
    dom2 != "NTr" ~ "str",
  )) %>% count(str_unstr_ck, str_unstr_ckr) 

data <- data %>%
  unite(col = type, c(str_unstr_ck, str_unstr_ckr), sep = "_")

data <- data %>% mutate(pct = n / sum(n))

data %>%
  ggplot(aes(x = "", y= pct)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

  
#############

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

# make matrix, choose upper triangle
# temp <- df
df <- df %>% pivot_wider(names_from = pdb1, values_from = rin_pid)
df <- as.matrix(df)
rownames(df) <- df[,1]
df <- df[, -1]

# remake data frame
df[upper.tri(df)]<-NA
df <- as.data.frame(df)
df$a <- rownames(df)
df <- df %>% pivot_longer(cols = 1:(ncol(df)-1), names_to = "b", values_to = "n"  )

# remove NAs
df <- df %>% drop_na(n)

# make numeric
df$n <- as.numeric(df$n)

# remove identical pairwise comparisons
df <- df %>% filter(a != b)

# pairwise RMSD & join
rmsd <- data.frame(a = c("7o7f", "7f1t", "5wb2", "6lfo", "6wwz", "7xa3", "7vl9"), b = c("5uiw", "5uiw", "4xt1", "5uiw", "5uiw", "5uiw", "5uiw"), rmsd = c(1.241, 0.905, 0.634, 1.889, 2.069, 1.529, 1.536))

df <- left_join(df, rmsd)
df %>%
  ggplot(aes(rmsd, 1)) +
  geom_jitter(shape = 21, colour = "black", fill = "white", size = 8, stroke = 0.5, width=0.) +
  # ylim(-0.05, 0.5) +
  # geom_hline(yintercept = c(0.1, 0.2)) +
  theme_minimal()



df <- df %>% filter(!(a %in% c("zheng", "7o7f", "7f1r", "4xt1", "5wb2", "4rws", "7sk3")))
df <- df %>% filter(!(b %in% c("zheng", "7o7f", "7f1r", "4xt1", "5wb2", "4rws", "7sk3")))



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
df$ck1 <- lookup$chemokine[match(unlist(df$a), lookup$pdbs)]
df$ckr1 <- lookup$receptor[match(unlist(df$a), lookup$pdbs)]
df$ck2 <- lookup$chemokine[match(unlist(df$b), lookup$pdbs)]
df$ckr2 <- lookup$receptor[match(unlist(df$b), lookup$pdbs)]
df <- df %>% unite(ck1_ck2, c(ck1, ck2))
df <- df %>% unite(ckr1_ckr2, c(ckr1, ckr2))


# plot
df$arb <- c(1)
df %>%
  # filter(a == "7sk3" | b == "7sk3") %>%
  ggplot(aes(arb, n)) +
  # geom_violin(trim = FALSE) +
  geom_jitter(shape = 21, colour = "black", fill = "white", size = 8, stroke = 0.5) +
  ylim(-0.05, 0.5) +
  geom_hline(yintercept = c(0.1, 0.2)) +
  theme_minimal()


#----
# Get CK and CKR no_pdb vs cons

# chemokine
data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  select(source_gnccn, dom1, no_contacts_all_ck, all_para_ck) %>%
  unique()

data <- data %>% mutate(str_unstr = case_when(
  dom1 == "NTc" ~ "unstr",
  dom1 != "NTc" ~ "str",
))

data %>%
  ggplot(aes(str_unstr, all_para_ck)) +
  geom_boxplot() +
  geom_violin() +
  theme_minimal()

data %>%
  filter(!(is.na(source_gnccn))) %>%
  ggplot(aes(str_unstr, no_contacts_all_ck)) +
  geom_boxplot() +
  # geom_violin() +
  theme_minimal()

data %>%
  filter(!(is.na(source_gnccn))) %>%
  # filter(dom1 == "NTc") %>%
  ggplot(aes(no_contacts_all_ck, all_para_ck, fill = str_unstr)) +
  geom_jitter(shape = 21, colour = "black", size = 8, stroke = 0.5) +
  # xlim(0,9) +
  ylim(0,1) +
  theme_minimal()

# receptor
data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  select(target_gnccn, dom2, no_contacts_all_ckr, all_para_ckr) %>%
  unique()


data <- data %>% mutate(str_unstr = case_when(
  dom2 == "NTr" ~ "unstr",
  dom2 != "NTr" ~ "str",
))


data %>%
  ggplot(aes(str_unstr, all_para_ckr)) +
  geom_boxplot() +
  geom_violin() +
  theme_minimal()

data %>%
  ggplot(aes(str_unstr, no_contacts_all_ckr)) +
  geom_boxplot() +
  # geom_violin() +
  theme_minimal()


data %>%
  filter(!(is.na(target_gnccn))) %>%
  # filter(dom2 == "NTr") %>%
  ggplot(aes(no_contacts_all_ckr, all_para_ckr, fill = str_unstr)) +
  geom_jitter(shape = 21, colour = "black", size = 8, stroke = 0.5) +
  xlim(0,9) +
  ylim(0,1) +
  theme_minimal()
