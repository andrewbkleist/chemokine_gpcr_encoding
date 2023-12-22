source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CHEMOKINE N-TERM ---------------------------------------------------------
total <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv")
frag <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
  select(motif) %>%
  unique()
slim <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
  filter(pct_ortho >= 0.5) %>%
  select(motif) %>%
  unique()
paste0("There are ", nrow(total), " total and ",
       nrow(frag), " unique fragments in chemokine N-termini,",
       " among which ", nrow(slim), "(", round(nrow(slim)/nrow(total),2), ")",
       " are putative SLiMs") 
rm(total, frag, slim)

# (2) GPCR N-TERM --------------------------------------------------------------
total <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
frag <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv") %>%
  select(motif) %>%
  unique()
slim <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv") %>%
  filter(pct_ortho >= 0.5) %>%
  select(motif) %>%
  unique()
paste0("There are ", nrow(total), " total and ",
       nrow(frag), " unique fragments in receptor N-termini,",
       " among which ", nrow(slim), "(", round(nrow(slim)/nrow(total),2), ")",
       " are putative SLiMs") 
rm(total, frag, slim)

# (3) GPCR N-TERM --------------------------------------------------------------
total <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv")
frag <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
  select(motif) %>%
  unique()
slim <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
  filter(pct_ortho >= 0.5) %>%
  select(motif) %>%
  unique()
paste0("There are ", nrow(total), " total and ",
       nrow(frag), " unique fragments in receptor ECL2,",
       " among which ", nrow(slim), "(", round(nrow(slim)/nrow(total),2), ")",
       " are putative SLiMs") 
rm(total, frag, slim)
