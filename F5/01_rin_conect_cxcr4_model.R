# Name:     02_rin_conect.R
# Updated:  20210516
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
  

##### 1: CXCL12 AND CXCR4 MOTIF CONTACTS #######################################
  
  # Define contacts of interest
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") 
  rin <- rin %>% filter(file == "ngo") %>%
    filter(source_gnccn %in% c("NTc.Cm8", "NTc.Cm7", "NTc.Cm6") | target_gnccn %in% c("NTr.Cm6", "NTr.Cm7", "NTr.Cm8")) %>%
    select(source_gnccn, target_gnccn) %>% unique()
  
  # (1.2) Write CONECT records for PDBs ----------------------------------------
  pdb.ngo <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
                          rin,
                          "ngo",
                          "01_structure_contacts/data/pdbs/ngo_model_clean.pdb",
                          "54_cxcr4_sat_mut/output/cxcr4_motif_rins.csv")
  
