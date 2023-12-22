source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import contacts
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(no_pdb == 1) %>%
  dplyr::select(source_gnccn, target_gnccn) %>% unique()

# (1) 5uiw
pdb.5uiw <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "5uiw",
                              "data/structure/raw/pdb_files/5uiw_ck_clean.pdb",
                              "output/F4S2/conect/5uiw_unique.csv")

# (2) 4rws
pdb.4rws <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "4rws",
                              "data/structure/raw/pdb_files/4rws_ck_clean.pdb",
                              "output/F4S2/conect/4rws_unique.csv")

# (3) 4xt1
pdb.4xt1 <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "4xt1",
                              "data/structure/raw/pdb_files/4xt1_ck_clean.pdb",
                              "output/F4S2/conect/4xt1_unique.csv")

# (4) 5wb2
pdb.5wb2 <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "5wb2",
                              "data/structure/raw/pdb_files/5wb2_clean.pdb",
                              "output/F4S2/conect/5wb2_unique.csv")

# (5) 6lfo
pdb.6lfo <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "6lfo",
                              "data/structure/raw/pdb_files/6lfo_clean.pdb",
                              "output/F4S2/conect/6lfo_unique.csv")


# (6) 6wwz
pdb.6wwz <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "6wwz",
                              "data/structure/raw/pdb_files/6wwz_clean.pdb",
                              "output/F4S2/conect/6wwz_unique.csv")

# (7) ngo
pdb.ngo <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                             rin,
                             "ngo",
                             "data/structure/raw/pdb_files/ngo_model_clean.pdb",
                             "output/F4S2/conect/ngo_unique.csv")

# (8) zheng
pdb.zheng <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                               rin,
                               "zheng",
                               "data/structure/raw/pdb_files/zheng_model_clean.pdb",
                               "output/F4S2/conect/zheng_unique.csv")

# (9) 7xbx
pdb.7xbx <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7xbx",
                              "data/structure/raw/pdb_files/7xbx_clean.pdb",
                              "output/F4S2/conect/7xbx_unique.csv")

# (10) 7vl9
pdb.7vl9 <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7vl9",
                              "data/structure/raw/pdb_files/7vl9_clean.pdb",
                              "output/F4S2/conect/7vl9_unique.csv")

# (11) 7sk3
pdb.7sk3 <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7sk3",
                              "data/structure/raw/pdb_files/7sk3_clean.pdb",
                              "output/F4S2/conect/7sk3_unique.csv")

# (12) 7o7f
pdb.7o7f <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7o7f",
                              "data/structure/raw/pdb_files/7o7f_clean.pdb",
                              "output/F4S2/conect/7o7f_unique.csv")

# (13) 7f1t
pdb.7f1t <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7f1t",
                              "data/structure/raw/pdb_files/7f1t_clean.pdb",
                              "output/F4S2/conect/7f1t_unique.csv")

# (14) 7f1r
pdb.7f1r <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7f1r",
                              "data/structure/raw/pdb_files/7f1r_clean.pdb",
                              "output/F4S2/conect/7f1r_unique.csv")

# (15) 7xa3
pdb.7xa3 <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "7xa3",
                              "data/structure/raw/pdb_files/7xa3_clean.pdb",
                              "output/F4S2/conect/7xa3_unique.csv")

# (16) 7xa3
pdb.8ic0 <- WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                              rin,
                              "8ic0",
                              "data/structure/raw/pdb_files/8ic0_clean.pdb",
                              "output/F4S2/conect/8ic0_unique.csv")

## LAST WRITTEN 20231117
# Outputs manually pasted onto ends of "cleaned" PDB files in same folder;
# Spaces manually adjusted to match CONECT record conventions