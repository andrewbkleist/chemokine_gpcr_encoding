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
    data <- GetContactList(i, j)
    temp <- data.frame(data, paste0(i), paste0(j))
    df <- rbind(df, temp)
  }
}
rm(temp, i, j, data)
colnames(df) <- c("source_gnccn", "target_gnccn", "type", "pdb1", "pdb2")

# remove self-by-self
df <- df %>% filter(pdb1 != pdb2)

# remove duplicates
temp <- df %>% dplyr::select(pdb1, pdb2)
temp <- data.frame(t(apply(temp,1,sort)))
temp <- temp[!duplicated(temp),]
colnames(temp) <- c("pdb1","pdb2")
temp$temp <- c("temp")
df <- left_join(df, temp)
df <- df %>% filter(temp == "temp") %>% dplyr::select(-temp)
rm(temp)

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
df$ck1 <- lookup$chemokine[match(unlist(df$pdb1), lookup$pdbs)]
df$ckr1 <- lookup$receptor[match(unlist(df$pdb1), lookup$pdbs)]
df$ck2 <- lookup$chemokine[match(unlist(df$pdb2), lookup$pdbs)]
df$ckr2 <- lookup$receptor[match(unlist(df$pdb2), lookup$pdbs)]
df <- df %>% unite(ck1_ck2, c(ck1, ck2))
df <- df %>% unite(ckr1_ckr2, c(ckr1, ckr2))
rm(lookup)

# add domain designations
lookup.domain <- read.csv("data/lookup/lookup_gnccn_to_domain.csv")
df$dom1 <- lookup.domain$dom[match(unlist(df[ ,"source_gnccn"]), lookup.domain$bwccn)]
df$dom2 <- lookup.domain$dom[match(unlist(df[ ,"target_gnccn"]), lookup.domain$bwccn)]  
rm(lookup.domain)

# reorder cols
df <- df %>% dplyr::select(source_gnccn, target_gnccn, dom1, dom2, type, 
                           pdb1, pdb2, ck1_ck2, ckr1_ckr2)

# write
# write_csv(df, "data/structure/processed/RIN_pairwise.csv")
# written 20231118