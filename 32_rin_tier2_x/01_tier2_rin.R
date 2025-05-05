# Name:     01_tier2_rin.R
# Updated:  20201104
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
  
##### 1: CC-/CXC- SPECIFIC #####################################################
  
  # (1.1) Define CC- and CXC-specific consensus networks -----------------------
    # import, clean, unite
    rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
      #filter(no_pdb >= 2) %>%
      select(file, source_gnccn, target_gnccn) 
    
    # similarities and differences among complexes
    pdb.5uiw <- rin %>% filter(file == "5uiw") %>% select(-file)
    pdb.6wwz <- rin %>% filter(file == "6wwz") %>% select(-file)
    pdb.4rws <- rin %>% filter(file == "4rws") %>% select(-file)
    pdb.6lfo <- rin %>% filter(file == "6lfo") %>% select(-file)
    pdb.ngo <- rin %>% filter(file == "ngo") %>% select(-file)
    pdb.4xt1 <- rin %>% filter(file == "4xt1") %>% select(-file)
    pdb.5wb2 <- rin %>% filter(file == "5wb2") %>% select(-file)
    pdb.zheng <- rin %>% filter(file == "zheng") %>% select(-file)
    
    
    # intersection CXC - CXCL12 model, CXCL8
    intersect.cxc <- intersect(pdb.6lfo, pdb.ngo)
    
    # intersection CC - CCL20, CCL5, Zheng model
    intersect.cc <- intersect(pdb.5uiw, pdb.6wwz)
    intersect.cc <- intersect(intersect.cc, pdb.zheng)
    
    # vMIP
    # intersect.cc.vmip  <- intersect(intersect.cc, pdb.4rws) # 8 (more CC-like)
    # intersect.cxc.vmip  <- intersect(intersect.cxc, pdb.4rws) # 4
    
    # setdiff
    setdiff.cc <- setdiff(intersect.cc, pdb.6lfo)
    setdiff.cc <- setdiff(setdiff.cc, pdb.ngo)
    
    setdiff.cxc <- setdiff(intersect.cxc, pdb.5uiw)
    setdiff.cxc <- setdiff(setdiff.cxc, pdb.6wwz)
    setdiff.cxc <- setdiff(setdiff.cxc, pdb.zheng)
    
    # intsersect
    intersect.cc.cxc <- intersect(intersect.cc, intersect.cxc)
    
    # remove
    rm(pdb.5uiw, pdb.6wwz, pdb.4rws, pdb.6lfo, pdb.ngo, pdb.4xt1, pdb.5wb2, pdb.zheng, rin)
  
  # (1.2) Write CONECT ---------------------------------------------------------
    # setdiff.cc <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
    #                         setdiff.cc,
    #                         "5uiw",
    #                         "01_structure_contacts/data/pdbs/5uiw_ck_clean.pdb",
    #                         "32_rin_tier2/output/5uiw_cc_specific.csv")
    
    # setdiff.cxc <- WriteCONECTcustom("01_structure_contacts/output/RIN_residue.csv",
    #                                 setdiff.cxc,
    #                                 "6lfo",
    #                                 "01_structure_contacts/data/pdbs/6lfo_clean.pdb",
    #                                 "32_rin_tier2/output/6lfo_cxc_specific.csv")
  
  # (1.3) Define consensus network (all structure "INTERSECTION") --------------
    
    # consensus
    # rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    #   filter(no_pdb >= 5) %>%
    #   select(source_gnccn, target_gnccn) %>% unique()
    # rin$type <- c("consensus")
    
    
  # (1.4) Write network for cytoscape ------------------------------------------
    
    # add label
    setdiff.cc$type <- c("cc")
    setdiff.cxc$type <- c("cxc")
    intersect.cc.cxc$type <- c("shared")
  
    # bind
    master <- bind_rows(setdiff.cc, setdiff.cxc, intersect.cc.cxc)  
    
    # write
    write_csv(master, "32_rin_tier2/output/cc_cxc_specific_cyto.csv")

       
    