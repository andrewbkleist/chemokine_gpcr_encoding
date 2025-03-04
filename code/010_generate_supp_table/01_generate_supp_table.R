source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) COMBINE, EXPORT GNOMAD
ck <- read_csv("data/variant/gnomad/processed/CK_GNOMAD_TABLE.csv") %>%
  select(-aln)
colnames(ck)[49] <- c("ccn_or_crn")

ckr <- read_csv("data/variant/gnomad/processed/CKR_GNOMAD_TABLE.csv") %>%
  select(-aln)
colnames(ckr)[49] <- c("ccn_or_crn")

data <- rbind(ck,ckr)
rm(ck,ckr)
data <- data %>% select(gene_symbol, ccn_or_crn,aa_ref,aa_pos,aa_alt,aln_pos, everything())

write_csv(data, "output/supp_tables/table_S3.csv")
rm(data)

# (2) COMBINE, EXPORT TCGA
ck <- read_csv("data/variant/tcga/processed/CK_TCGA_COUNTS.csv")
colnames(ck)[2] <- c("ccn_or_cgn")
ckr <- read_csv("data/variant/tcga/processed/CKR_TCGA_COUNTS.csv")
colnames(ckr)[2] <- c("ccn_or_cgn")
ckr <- ckr %>% separate(col = ccn_or_cgn, into = c("a", "gn"), sep = "gn", remove = FALSE) %>%
  select(protein,gn,cancer_mut_count)
colnames(ckr)[2] <- c("ccn_or_cgn")

data <- rbind(ck,ckr)
rm(ck,ckr)

write_csv(data, "output/supp_tables/table_S4.csv")
rm(data)


