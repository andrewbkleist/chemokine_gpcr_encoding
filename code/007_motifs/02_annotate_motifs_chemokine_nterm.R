# Imports motif list and adds annotations (e.g. count among orthologs, paralogs,
# paralog conservation, etc)

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) FILTERING & CLEANING -----------------------------------------------------
data <- read_csv("data/motif/processed/CK_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")

# REMOVE DOUBLES
data <- data %>% dplyr::select(-motif_no) %>% distinct() 

# ADD NUMBER OF ORTHOLOGS
no.ortho <- read.csv("data/lookup/no_ortho.csv")
data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$ck)]
data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$ck)]
rm(no.ortho)

# REPLACE "NA" WITH ASN-ALA
data$motif[is.na(data$motif)] <- c("AsnAla")

# REMOVE XCL1, XCL2, CCL4L1
data <- data %>% filter(protein != "xcl1" & protein != "xcl2" & protein != "ccl4l1")

# (2) CONSERVATION SCORING AT DIFFERENT LEVELS (ORTHOLOGS, FAMILY, ETC) --------
# INDIVIDUAL - motif occurance across a single chemokine
indiv <- data %>% dplyr::count(motif, protein, class, mer, mask, no_ortho) %>% 
  mutate(pct_ortho = n / no_ortho)
colnames(indiv)[6] <- c("total_ortho")
colnames(indiv)[7] <- c("count_ortho")
indiv <- indiv %>% dplyr::select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)

# DEFINE HUMAN SEQUENCES (counts will be done among human paralogs)
human <- data %>% filter(grepl("HUMAN", data$file))

# SUBFAMILY - motif occurrence across each class
# Motids are counted regardless of ortholog conservation. They are counted 
# according to the number of times they are represented in chemokine subfamilies 
# (ie CC, CXC) or the chemokine family (ie all chemokines).Evaluation occurs on 
# paralog level

family <- human %>% dplyr::count(motif, class)
colnames(family)[3] <- c("count_family")

freq <- left_join(indiv, family)
freq[is.na(freq)] <- 0

# SUPERFAMILY - motif occurance across all chemokines
super <- human %>% dplyr::count(motif)
colnames(super)[2] <- c("count_super")

freq <- left_join(freq, super)
freq[is.na(freq)] <- 0

# ADD PCT COLUMNS
freq <- freq %>%
  mutate(total_family = case_when(
    class == "cc" ~ 27, # recall that CCL4L1 has been removed, so 27 instead of 28
    class == "cxc" ~ 17,
    class == "cx3c" ~ 1 # note that XC family has been removed
  ))

freq <- freq %>% mutate(pct_family = count_family / total_family)
freq$total_super <- 43
freq <- freq %>% mutate(pct_super = count_super / total_super)

# (3) CLEAN AND WRITE ----------------------------------------------------------
freq <- freq %>% dplyr::select(motif, protein, class, mer, mask, 
                        count_ortho, total_ortho, pct_ortho,
                        count_family, total_family, pct_family,
                        count_super, total_super, pct_super)

# add column for whether found in human or not, filter for human only motifs
human <- human %>% dplyr::select(protein, motif)
human$human <- c(1)
freq <- left_join(freq, human)
freq <- freq %>% filter(human == 1)
freq <- freq %>% dplyr::select(-human)

# write csv
# write_csv(freq, "data/motif/processed/CK_MOTIF_FREQUENCY.csv")
# WRITTEN 20231001
rm(data, family, freq, indiv, super, human)