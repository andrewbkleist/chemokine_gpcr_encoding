source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# define all PDBs
pdbs <- c("5uiw", "7o7f", "7f1r", "zheng",
          "7vl9", "7xa3", "7f1t", "6wwz",
          "8ic0", "6lfo", "ngo",
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


# add chemokine, receptor information
chemokine <- c("CCL5", "CCL5", "CCL5", "CCL5", 
               "CCL15", "CCL2", "CCL3", "CCL20", 
               "CXCL8", "CXCL8", "CXCL12", "CXCL12",
               "CX3CL1", "vMIPII", "CX3CL1", "CX3CL1"
)

receptor <- c("CCR5", "CCR5", "CCR5", "CCR5", 
              "CCR1", "CCR2", "CCR5", "CCR6",
              "CXCR1", "CXCR2", "CXCR4",
              "ACKR3", "CX3CR1",
              "CXCR4", "US28", "US28")
lookup <- data.frame(pdbs, chemokine, receptor)

# map chemokine/receptor names to df, unite columns
df$ck1 <- lookup$chemokine[match(unlist(df$a), lookup$pdbs)]
df$ckr1 <- lookup$receptor[match(unlist(df$a), lookup$pdbs)]
df$ck2 <- lookup$chemokine[match(unlist(df$b), lookup$pdbs)]
df$ckr2 <- lookup$receptor[match(unlist(df$b), lookup$pdbs)]

# write csv
write_csv(df, "data/structure/processed/RIN_all_by_all_pid.csv")
