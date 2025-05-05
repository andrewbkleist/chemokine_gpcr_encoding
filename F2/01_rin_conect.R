# Name:     01_rin_conect.R
# Updated:  20221007
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)


##### FUNCTIONS ################################################################

  WriteCONECTcustom <- function(RINFILE, RINDF, PDBID, PDBFILE, OUTPUT){
    
    # import rinfile
    rin <- read_csv(RINFILE) %>% filter(file == PDBID)
    rin.df <- RINDF
    rin.df$sele <- c("yes")
    rin <- left_join(rin, rin.df)
    rin <- rin %>% filter(sele == "yes")
    rin <- rin %>% select(-sele)
    
    # read PDB, make df, select relevant columns
    pdb <- read.pdb(PDBFILE)
    pdb_df <- as.data.frame(pdb$atom)
    pdb_conv <- pdb_df %>% select(chain, resno, elety, eleno)
    pdb_conv <- pdb_conv %>% filter(elety == "CA")
    ck <- pdb_conv %>% filter(chain == "A")
    ckr <- pdb_conv %>% filter(chain == "B")
    
    # map atom indices to RIN file
    rin$ca1 <- ck$eleno[match(unlist(rin$ResNum1), ck$resno)]
    rin$ca2 <- ckr$eleno[match(unlist(rin$ResNum2), ckr$resno)]
    
    # clean up and write
    rin$CONECT <- c("CONECT")
    rin <- rin %>% select(CONECT, ca1, ca2)
    write_csv(rin, OUTPUT)
    
    # return
    return(rin)
    
    # remove
    rm(rin, pdb, pdb_df, ck, ckr, pdb_conv)
  }
  
  
  
##### 1: ALL CONTACTS, ALL STRUCTURES ##########################################
  
  # import contacts
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    # filter(no_pdb >= 5) %>%
    select(source_gnccn, target_gnccn) %>% unique()
  
  # # (1) 5uiw
  # pdb.5uiw <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                         rin,
  #                         "5uiw",
  #                         "01_structure_contacts/data/pdbs/5uiw_ck_clean.pdb",
  #                         "20_rin_alone/output/5uiw_all_rins.csv")
  # 
  # # (2) 4rws
  # pdb.4rws <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "4rws",
  #                               "01_structure_contacts/data/pdbs/4rws_ck_clean.pdb",
  #                               "20_rin_alone/output/4rws_all_rins.csv")
  # 
  # # (3) 4xt1
  # pdb.4xt1 <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "4xt1",
  #                               "01_structure_contacts/data/pdbs/4xt1_ck_clean.pdb",
  #                               "20_rin_alone/output/4xt1_all_rins.csv")
  # 
  # # (4) 5wb2
  # pdb.5wb2 <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "5wb2",
  #                               "01_structure_contacts/data/pdbs/5wb2_clean.pdb",
  #                               "20_rin_alone/output/5wb2_all_rins.csv")
  # 
  # # (5) 6lfo
  # pdb.6lfo <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "6lfo",
  #                               "01_structure_contacts/data/pdbs/6lfo_clean.pdb",
  #                               "20_rin_alone/output/6lfo_all_rins.csv")
  # 
  # 
  # # (6) 6wwz
  # pdb.6wwz <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "6wwz",
  #                               "01_structure_contacts/data/pdbs/6wwz_clean.pdb",
  #                               "20_rin_alone/output/6wwz_all_rins.csv")
  # 
  # # (7) ngo
  # pdb.ngo <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                              rin,
  #                              "ngo",
  #                              "01_structure_contacts/data/pdbs/ngo_model_clean.pdb",
  #                              "20_rin_alone/output/ngo_all_rins.csv")
  # 
  # # (8) zheng
  # pdb.zheng <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                                rin,
  #                                "zheng",
  #                                "01_structure_contacts/data/pdbs/zheng_model_clean.pdb",
  #                                "20_rin_alone/output/zheng_all_rins.csv")
  # 
  # # (9) 7xbx
  # pdb.7xbx <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                                rin,
  #                                "7xbx",
  #                                "01_structure_contacts/data/pdbs/7xbx_clean.pdb",
  #                                "20_rin_alone/output/7xbx_all_rins.csv")
  # 
  # # (10) 7vl9
  # pdb.7vl9 <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "7vl9",
  #                               "01_structure_contacts/data/pdbs/7vl9_clean.pdb",
  #                               "20_rin_alone/output/7vl9_all_rins.csv")
  # 
  # # (11) 7sk3
  # pdb.7sk3 <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "7sk3",
  #                               "01_structure_contacts/data/pdbs/7sk3_clean.pdb",
  #                               "20_rin_alone/output/7sk3_all_rins.csv")
  # 
  # # (12) 7o7f
  # pdb.7o7f <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "7o7f",
  #                               "01_structure_contacts/data/pdbs/7o7f_clean.pdb",
  #                               "20_rin_alone/output/7o7f_all_rins.csv")
  # 
  # # (13) 7f1t
  # pdb.7f1t <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "7f1t",
  #                               "01_structure_contacts/data/pdbs/7f1t_clean.pdb",
  #                               "20_rin_alone/output/7f1t_all_rins.csv")
  # 
  # # (14) 7f1r
  # pdb.7f1r <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "7f1r",
  #                               "01_structure_contacts/data/pdbs/7f1r_clean.pdb",
  #                               "20_rin_alone/output/7f1r_all_rins.csv")
  # 
  # # (15) 7xa3
  # pdb.7xa3 <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                               rin,
  #                               "7xa3",
  #                               "01_structure_contacts/data/pdbs/7xa3_clean.pdb",
  #                               "20_rin_alone/output/7xa3_all_rins.csv")
  
  ## LAST WRITTEN 20221007
  

##### 2: CONSENSUS NETWORK (>5 PDBs) ###########################################
  
  # (1.1) Define consensus network (all structure "INTERSECTION") --------------
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(no_pdb >= 5) %>%
    select(source_gnccn, target_gnccn) %>% unique()

  # (1.2) Write CONECT records for PDBs ----------------------------------------
  # pdb.5uiw <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
  #                         rin,
  #                         "5uiw", 
  #                         "01_structure_contacts/data/pdbs/5uiw_ck_clean.pdb", 
  #                         "30_rin_analysis/output/5uiw_consensus.csv")
  
  # (1.3) Write contacts for later analysis ------------------------------------
  rin$type <- c("consensus")
  write_csv(rin, "20_rin_alone/output/rin_consensus_5uiw.csv")
  
