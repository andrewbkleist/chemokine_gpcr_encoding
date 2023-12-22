source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import contacts
rin.5uiw <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file == "8ic0") %>%
  dplyr::select(source_gnccn, target_gnccn) %>% unique()

rin.7o7f <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file == "6lfo") %>%
  dplyr::select(source_gnccn, target_gnccn) %>% unique()

rin.5uiw.not.7o7f <- dplyr::setdiff(rin.5uiw, rin.7o7f)
rin.7o7f.not.5uiw <- dplyr::setdiff(rin.7o7f, rin.5uiw)
rin.7o7f.and.5uiw <- dplyr::intersect(rin.7o7f, rin.5uiw)

paste0("5uiw has ", nrow(rin.5uiw.not.7o7f), " unique contacts,",
       " 7o7f has ", nrow(rin.7o7f.not.5uiw), " unique contacts,",
       " and the files share ", nrow(rin.7o7f.and.5uiw), " unique contacts")

# DID NOT RE-WRITE CONECT FILE, WRITTEN PREVIOUSLY IN A DIFFERENT FOLDER DIRECTORY
# BUT NO CHANGES SINCE THEN. FOR ILLUSTRATION ONLY
# # (1) 5uiw not 7o7f
# rin.5uiw.not.7o7f <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
#                               rin.5uiw.not.7o7f,
#                               "5uiw",
#                               "01_structure_contacts/data/pdbs/5uiw_ck_clean.pdb",
#                               "50_rin_alone/output/5uiw_not_7o7f_rins.csv")
# 
# # (2) 7o7f not 5uiw
# rin.7o7f.not.5uiw <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
#                               rin.7o7f.not.5uiw,
#                               "7o7f",
#                               "01_structure_contacts/data/pdbs/7o7f_clean.pdb",
#                               "50_rin_alone/output/7o7f_not_5uiw_rins.csv")
