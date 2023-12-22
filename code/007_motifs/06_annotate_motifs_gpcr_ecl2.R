source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) FILTERING & CLEANING -----------------------------------------------------
# setup data
data1 <- read_csv("data/motif/processed/CKR_UNSTRUCTURED_ECL2A_ALL_MERS_CYSLESS.csv")
data2 <- read_csv("data/motif/processed/CKR_UNSTRUCTURED_ECL2B_ALL_MERS_CYSLESS.csv")
data <- bind_rows(data1, data2)
rm(data1, data2)
data <- data %>% dplyr::select(-motif_no) %>% distinct() 
# remove doubles of motifs within same sequence

no.ortho <- read.csv("data/lookup/no_ortho.csv")
no.ortho$protein <- as.character(no.ortho$protein)
data$no_ortho <- no.ortho$no_ortho[match(unlist(data$protein), no.ortho$protein)]
data$partners <- no.ortho$partners[match(unlist(data$protein), no.ortho$protein)]
rm(no.ortho)
data$motif[is.na(data$motif)] <- c("AsnAla") # no gaps are included so all NAs must be AsnAla

# (2) CONSERVATION SCORING AT DIFFERENT LEVELS (ORTHOLOGS, FAMILY, ETC) --------
# INDIVIDUAL - motif occurance across a single chemokine
indiv <- data %>% dplyr::count(motif, protein, class, mer, mask, no_ortho) %>% 
  mutate(pct_ortho = n / no_ortho)
colnames(indiv)[6] <- c("total_ortho")
colnames(indiv)[7] <- c("count_ortho")
indiv <- indiv %>% dplyr::select(motif, protein, class, mer, mask, count_ortho, total_ortho, pct_ortho)

# Note that NA --> AsnAla and doubles already addressed

# DEFINE HUMAN SEQUENCES (counts will be done among human paralogs)
human <- data %>% filter(grepl("human", data$file))

# SUBFAMILY - motif occurrence across each class
# This is changed since original submission; motifs are counted regardless
# of ortholog conservation. Instead they are counted according to the
# number of times they are represented in receptor subfamilies (ie CC, CXC)
# or the receptor family (ie all chemokine receptors).
# Evaluation occurs on paralog level, not abundance of ortholog sequences
# due to assymetric data sets for each receptor (ie some have 50 sequences)
family <- human %>% dplyr::count(motif, class)
colnames(family)[3] <- c("count_family")

freq <- left_join(indiv, family)
freq[is.na(freq)] <- 0


# FAMILY - motif occurance across all receptors
super <- human %>% dplyr::count(motif)
colnames(super)[2] <- c("count_super")

freq <- left_join(freq, super)
freq[is.na(freq)] <- 0

# ADD PCT COLUMNS
freq <- freq %>%
  mutate(total_family = case_when(
    class == "cc" ~ 10,
    class == "cxc" ~ 6,
    class == "ack" ~ 5,
    class == "xc" ~ 1,
    class == "cx3c" ~ 1
  ))
freq <- freq %>% mutate(pct_family = count_family / total_family)
freq$total_super <- 23
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

# write output
# write_csv(freq, "data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv")
# WRITTEN 20231127
rm(data, family, freq, indiv, super, human)
