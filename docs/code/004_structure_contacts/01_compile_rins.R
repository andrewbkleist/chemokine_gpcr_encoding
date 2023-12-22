# Imports Protein Contact Atlas output files with residue interaction networks
# (RINs), adds common numbering (CCN and CGN) and SSE information, and compiles 
# into *.csv file as output. Retains intra- and intermolecular RINs, reduces to
# one RIN per residue pair

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

  # load CCN and BW conversion file
  lookup.bwccn <- read.csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv")
  
  # (1) CXCL8:CXCR1
  cna.1ilp <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/1ILPCKCLEAN1538673510.txt",
                        "clean_1ilp_ck", "ccn_1ilp_ck",     # object resno, ccn
                        "clean_1ilp_ckr", "bw_1ilp_ckr",    # ligand resno, bw
                        "1ilp",                             # PDB ID
                        "soluble")                          # PDB "class"

  
  # (2) CCL11:CCR3
  cna.2mpm <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/2MPMCKCLEAN1538684712.txt",
                          "clean_2mpm_ck", "ccn_2mpm_ck",   # object resno, ccn
                          "clean_2mpm_ckr", "bw_2mpm_ckr",  # ligand resno, bw
                          "2mpm",                           # PDB ID
                          "nmr")                            # PDB "class"
  
  # (3) CXCL12:CXCR4
  cna.2n55 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/2N55CKCLEAN1538684756.txt",
                          "clean_2n55_ck", "ccn_2n55_ck",     # object resno, ccn
                          "clean_2n55_ckr", "bw_2n55_ckr",  # ligand resno, bw
                          "2n55",                           # PDB ID
                          "soluble")                        # PDB "class"
  
  # (4) vMIPII:CXCR4
  cna.4rws <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/4RWSCKCLEAN1538684823.txt",
                          "clean_4rws_ck", "ccn_4rws_ck",     # object resno, ccn
                          "clean_4rws_ckr", "bw_4rws_ckr",  # ligand resno, bw
                          "4rws",                           # PDB ID
                          "full")                           # PDB "class"
  
  # (5) CX3CL1:US28
  cna.4xt1 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/4XT1CKCLEAN1538684872.txt",
                          "clean_4xt1_ck", "ccn_4xt1_ck",     # object resno, ccn
                          "clean_4xt1_ckr", "bw_4xt1_ckr",  # ligand resno, bw
                          "4xt1",                           # PDB ID
                          "full")                           # PDB "class"
  
  # (6) CCL5[5P7]:CCR5
  cna.5uiw <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/5UIWCKCLEAN1538684915.txt",
                          "clean_5uiw_ck", "ccn_5uiw_ck",     # object resno, ccn
                          "clean_5uiw_ckr", "bw_5uiw_ckr",  # ligand resno, bw
                          "5uiw",                           # PDB ID
                          "full")                           # PDB "class"

  # (7) CXCL12:CXCR4
  cna.2k05 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/2K05CLEAN1538684566.txt",
                          "clean_2k05_ck", "ccn_2k05_ck",   # object resno, ccn
                          "clean_2k05_ckr", "bw_2k05_ckr",  # ligand resno, bw
                          "2k05",                           # PDB ID
                          "soluble")                        # PDB "class"
  
  # (8) CCL5:CCR5
  cna.6fgp <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/6FGPCLEAN1538688371.txt",
                        "clean_6fgp_ck", "ccn_6fgp_ck",     # object resno, ccn
                        "clean_6fgp_ckr", "bw_6fgp_ckr",    # ligand resno, bw
                        "6fgp",                             # PDB ID
                        "soluble")                          # PDB "class"
  
  # (9) CX3CL1.35:US28
  cna.5wb2 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/5WB2CLEAN1538684959.txt",
                        "clean_5wb2_ck", "ccn_5wb2_ck",     # object resno, ccn
                        "clean_5wb2_ckr", "bw_5wb2_ckr",    # ligand resno, bw
                        "5wb2",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (10) CCL20:CCR6
  cna.6wwz <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/6WWZCLEAN1598531253.txt",
                        "clean_6wwz_ck", "ccn_6wwz_ck",     # object resno, ccn
                        "clean_6wwz_ckr", "bw_6wwz_ckr",    # ligand resno, bw
                        "6wwz",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (11) CXCL8:CXCR2
  cna.6lfo <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/6LFOCLEAN1598542701.txt",
                        "clean_6lfo_ck", "ccn_6lfo_ck",     # object resno, ccn
                        "clean_6lfo_ckr", "bw_6lfo_ckr",    # ligand resno, bw
                        "6lfo",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (12) CXCL12:CXCR4 **MODEL**
  cna.ngo <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/NGOMODELCLEAN1600451722.txt",
                       "clean_ngo_ck", "ccn_ngo_ck",     # object resno, ccn
                       "clean_ngo_ckr", "bw_ngo_ckr",    # ligand resno, bw
                       "ngo",                            # PDB ID
                       "full")                           # PDB "class"
  
  # (13) CCL5:CCR5 **MODEL**
  cna.zheng <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/ZHENGMODELCLEAN1603803275.txt",
                         "clean_zheng_ck", "ccn_zheng_ck",     # object resno, ccn
                         "clean_zheng_ckr", "bw_zheng_ckr",    # ligand resno, bw
                         "zheng",                              # PDB ID
                         "full")                              # PDB "class"
  
  # (14) CX3CL1:CX3CR1
  cna.7xbx <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7XBXCLEAN1659321025.txt",
                        "clean_7xbx_ck", "ccn_7xbx_ck",     # object resno, ccn
                        "clean_7xbx_ckr", "bw_7xbx_ckr",    # ligand resno, bw
                        "7xbx",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (15) CCL15:CCR1
  cna.7vl9 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7VL9CLEAN1659321185.txt",
                        "clean_7vl9_ck", "ccn_7vl9_ck",     # object resno, ccn
                        "clean_7vl9_ckr", "bw_7vl9_ckr",    # ligand resno, bw
                        "7vl9",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (16) CXCL12:ACKR3
  cna.7sk3 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7SK3CLEAN1659321603.txt",
                        "clean_7sk3_ck", "ccn_7sk3_ck",     # object resno, ccn
                        "clean_7sk3_ckr", "bw_7sk3_ckr",    # ligand resno, bw
                        "7sk3",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (17) CCL5[6P4]:CCR5
  cna.7o7f <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7O7FCLEAN1659321637.txt",
                        "clean_7o7f_ck", "ccn_7o7f_ck",     # object resno, ccn
                        "clean_7o7f_ckr", "bw_7o7f_ckr",    # ligand resno, bw
                        "7o7f",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (18) CCL3:CCR5
  cna.7f1t <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7F1TCLEAN1659321674.txt",
                        "clean_7f1t_ck", "ccn_7f1t_ck",     # object resno, ccn
                        "clean_7f1t_ckr", "bw_7f1t_ckr",    # ligand resno, bw
                        "7f1t",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (19) CCL5:CCR5
  cna.7f1r <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7F1RCLEAN1659321718.txt",
                        "clean_7f1r_ck", "ccn_7f1r_ck",     # object resno, ccn
                        "clean_7f1r_ckr", "bw_7f1r_ckr",    # ligand resno, bw
                        "7f1r",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (20) gp120:CCR5
  cna.6meo <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/6MEOCLEAN1659334429.txt",
                        "clean_6meo_ck", "ccn_6meo_ck",     # object resno, ccn
                        "clean_6meo_ckr", "bw_6meo_ckr",    # ligand resno, bw
                        "6meo",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (21) CCL2:CCR2
  cna.7xa3 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/7XA3CLEAN1663721528.txt",
                        "clean_7xa3_ck", "ccn_7xa3_ck",     # object resno, ccn
                        "clean_7xa3_ckr", "bw_7xa3_ckr",    # ligand resno, bw
                        "7xa3",                             # PDB ID
                        "full")                             # PDB "class"
  
  
  # (22) CXCL8:CXCR1
  cna.8ic0 <- LoadAtlas("data/structure/raw/protein_contact_atlas_files/8IC0CLEAN1695496384.txt",
                        "clean_8ic0_ck", "ccn_8ic0_ck",     # object resno, ccn
                        "clean_8ic0_ckr", "bw_8ic0_ckr",    # ligand resno, bw
                        "8ic0",                             # PDB ID
                        "full")                             # PDB "class"
  
  # bind all contacts into a single df (9 PDBs), then remove
  cna.master <- bind_rows(cna.1ilp, cna.2mpm, cna.2n55, cna.4rws, cna.4xt1, 
                      cna.5uiw, cna.2k05, cna.6fgp, cna.5wb2, cna.6wwz, 
                      cna.6lfo, cna.ngo, cna.zheng, cna.7xbx, cna.7vl9, 
                      cna.7sk3, cna.7o7f, cna.7f1t, cna.7f1r, cna.6meo, 
                      cna.7xa3, cna.8ic0,
                      )

  rm(cna.1ilp, cna.2mpm, cna.2n55, cna.4rws, cna.4xt1, 
     cna.5uiw, cna.2k05, cna.6fgp, cna.5wb2, cna.6wwz, 
     cna.6lfo, cna.ngo, cna.zheng, cna.7xbx, cna.7vl9, 
     cna.7sk3, cna.7o7f, cna.7f1t, cna.7f1r, cna.6meo, 
     cna.7xa3, cna.8ic0)
  
  # add domain designations
  lookup.domain <- read.csv("data/lookup/lookup_gnccn_to_domain.csv")
  cna.master$dom1 <- lookup.domain$dom[match(unlist(cna.master[ ,"source_gnccn"]), lookup.domain$bwccn)]
  cna.master$dom2 <- lookup.domain$dom[match(unlist(cna.master[ ,"target_gnccn"]), lookup.domain$bwccn)]  
  rm(lookup.bwccn, lookup.domain)
  
  # remove unecessary columns
  cna.master <- cna.master %>% select(-PDB, -SS1, -SS2)
  
  # make unique such that each residue can only make one unique contact
  cna.master <- cna.master %>%
    select(-Number.of.atomic.contacts, -Atoms, -Chain.Types, -Distance) %>%
    unique()

  # write df
  # write_csv(cna.master, "data/structure/processed/RIN_residue.csv") # LAST WRITTEN 20230924
  
  