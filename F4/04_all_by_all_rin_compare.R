# Name:     04_all_by_all_rin_compare.R
# Updated:  20230501
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggridges)

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

# # remove identical pairwise comparisons
# # choose 1 CCL5:CCR5 example, remove viral-containing
df <- df %>% filter(a != b)
# df <- df %>% filter(!(a %in% c("zheng", "7o7f", "7f1r", "4xt1", "5wb2", "4rws", "7sk3")))
# df <- df %>% filter(!(b %in% c("zheng", "7o7f", "7f1r", "4xt1", "5wb2", "4rws", "7sk3")))
# 


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
rm(lookup)

# 
ccl5_rltd <- data.frame(pdb = c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", "7xa3", "7f1t","5uiw", "7o7f", "7f1r", "zheng", "7vl9", "7xa3", "7f1t"), type = "ccl5_rltd")
other_cc <- data.frame(pdb = c("6wwz"), type = "other_cc")
cxc <- data.frame(pdb = c("6lfo", "ngo"), type = "cxc")
other <- data.frame(pdb = c("7sk3", "7xbx"), type = "other")
viral <- data.frame(pdb = c("4rws", "4xt1", "5wb2"), type = "viral")
lookup <- rbind(ccl5_rltd,other_cc, cxc, other, viral)
rm(ccl5_rltd,other_cc, cxc, other, viral)

df$type1 <- lookup$type[match(unlist(df$a), lookup$pdb)]
df$type2 <- lookup$type[match(unlist(df$b), lookup$pdb)]
df <- df %>% unite(type12, c(type1, type2))
df <- df %>% 
  mutate(type12summ = case_when(
    type12 %in% c("ccl5_rltd_ccl5_rltd") ~ "intra_network_same_sub",
    type12 %in% c("cxc_cxc", "other_cc_ccl5_rltd") ~ "inter_network_same_sub",
    type12 %in% c("cxc_ccl5_rltd", "other_other_cc", "cxc_other_cc", "other_ccl5_rltd", "other_cxc", "other_other") ~ "diff_sub",
    type12 %in% c("viral_ccl5_rltd", "viral_other_cc", "viral_cxc", "viral_other", "viral_viral") ~ "viral"
  ))

# # add comparison-type column
# cc <- data.frame(pdb = c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", "7xa3", "7f1t", "6wwz"), type = "cc")
# cxc <- data.frame(pdb = c("6lfo", "ngo"), type = "cxc")
# other <- data.frame(pdb = c("7sk3", "7xbx"), type = "other")
# viral <- data.frame(pdb = c("4rws", "4xt1", "5wb2"), type = "viral")
# lookup <- rbind(cc, cxc, other, viral)
# rm(cc, cxc, other, viral)
# 
# df$type1 <- lookup$type[match(unlist(df$a), lookup$pdb)]
# df$type2 <- lookup$type[match(unlist(df$b), lookup$pdb)]
# df <- df %>% unite(type12, c(type1, type2))
# df <- df %>% 
#   mutate(type12summ = case_when(
#     type12 %in% c("cc_cc", "cxc_cxc") ~ "intra_sub",
#     type12 %in% c("cxc_cc", "other_cc", "other_cxc", "other_other") ~ "inter_sub",
#     type12 %in% c("viral_cc", "viral_cxc", "viral_other", "viral_viral") ~ "viral"
# ))


median(df$n)

# plot
df$arb <- c(1)
df %>%
  filter(type12summ != "viral") %>%
  filter(a == "ngo" & b == "7vl9") %>%
  ggplot(aes(arb, n, fill = type12summ)) +
  # geom_violin(trim = FALSE) +
  geom_point() +
  # geom_jitter(shape = 21, colour = "black",  size = 6, stroke = 0.5) +
  ylim(-0.05, 0.52) +
  geom_hline(yintercept = c(0.1, 0.2)) +
  scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  theme_minimal()

ggsave(filename = "rin_pct_compare_dot.pdf", 
       plot = last_plot(), path = "F4/output/",
       width = 3.5,
       height = 6)

df %>%
  filter(type12summ != "viral") %>%
  ggplot(aes(x = n, y = type12summ, fill = type12summ)) +
  geom_density_ridges() +
  scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  theme_minimal() 

ggsave(filename = "rin_pct_compare_ridge.pdf", 
       plot = last_plot(), path = "F4/output/",
       width = 8,
       height = 3)
