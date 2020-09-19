# Name:     01_cna_pca.r
# Updated:  20200918
# Author:   Andrew Kleist

# packages, working directory
library(here) # sets wd
library(tidyverse)

# NOTE 1 (20200918):
# all complexes were prepared by (1) removing hydrogens (via Bio3d clean.pdb);
# (2) the chemokine was made "Chain A" and the receptor "Chain B";
# (3) The "cleaned" PDB was submitted to Protein Contact Atlas 
# (https:/www.mrc-lmb.cam.ac.uk/rajini/index.html) with default settings
# GPCRdb nomenclature (GN) and Common Chemokine Numbering (CCN) designations
# were mapped using manually compiled "lookup table"; GN positions were modified
# to be inclusive of the conserved N-terminal Cys (designated 1x22);
# "Normal" GN positions from GPCRdb are identical past position 1x30

# NOTE 2 (20200918):
# Contact mappings were validated on 20200918 by spot checks...
# Chose 3 pairwise contacts - made sure spatially made sense & GN/CCN match
#     - 1ilp / 2k05 / 2mpm / 2n55 / 6fgp
#     - 4rws / 5uiw / 5wb2 / 6lfo / 6wwz / ngo model
#
# The following source/target/dom1/dom2 fields have NAs: 
# (1) 1ilp NH2 and ACA "atoms"
# (2) "GS" at the beginning of CXCR4 from 2n55

##### FUNCTIONS ################################################################

  # FUNCTION LoadAtlas: 
  
  LoadAtlasInter <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                          LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    # load file
    pca.object <- read.table(PCA.OUTPUT, sep="\t", header=TRUE)
    pca.object <- pca.object %>% filter(Chain1 == "A" & Chain2 =="B")
    
    # substitute for BW and/or CCN names
    object.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.OBJ.NAME][match(pca.object[ ,"ResNum1"], lookup.bwccn[ ,LOOK.OBJ.RES])])
    colnames(object.comm.name) <- c("source_gnccn")
    ligand.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.LIG.NAME][match(pca.object[ ,"ResNum2"], lookup.bwccn[ , LOOK.LIG.RES])])
    colnames(ligand.comm.name) <- c("target_gnccn")
    pca.object <- cbind(pca.object, object.comm.name)
    pca.object <- cbind(pca.object, ligand.comm.name)
    
    # add file name and PDB "class" (eg ck_ckr, cc_dimer, cxc_dimer)
    pdb.des <- as.data.frame(rep(PDBOBJ, nrow(pca.object)))
    colnames(pdb.des) <- c("file")
    pca.object <- cbind(pca.object, pdb.des)
    class.des <- as.data.frame(rep(CLASS, nrow(pca.object)))
    colnames(class.des) <- c("class")
    pca.object <- cbind(pca.object, class.des)
    
    # return final df, remove other objects
    return(pca.object)
    rm(pca.object, object.comm.name, ligand.comm.name, pca.object, pdb.des, class.des)
  }
  
  LoadAtlasCK <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                             LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    # load file
    pca.object <- read.table(PCA.OUTPUT, sep="\t", header=TRUE)
    pca.object <- pca.object %>% filter(Chain1 == "A" & Chain2 =="A")
    
    # substitute for BW and/or CCN names
    object.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.OBJ.NAME][match(pca.object[ ,"ResNum1"], lookup.bwccn[ ,LOOK.OBJ.RES])])
    colnames(object.comm.name) <- c("source_gnccn")
    ligand.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.OBJ.NAME][match(pca.object[ ,"ResNum2"], lookup.bwccn[ , LOOK.OBJ.RES])])
    colnames(ligand.comm.name) <- c("target_gnccn")
    pca.object <- cbind(pca.object, object.comm.name)
    pca.object <- cbind(pca.object, ligand.comm.name)
    
    # add file name and PDB "class" (eg ck_ckr, cc_dimer, cxc_dimer)
    pdb.des <- as.data.frame(rep(PDBOBJ, nrow(pca.object)))
    colnames(pdb.des) <- c("file")
    pca.object <- cbind(pca.object, pdb.des)
    class.des <- as.data.frame(rep(CLASS, nrow(pca.object)))
    colnames(class.des) <- c("class")
    pca.object <- cbind(pca.object, class.des)
    
    # return final df, remove other objects
    return(pca.object)
    rm(pca.object, object.comm.name, ligand.comm.name, pca.object, pdb.des, class.des)
  }
  
  
  LoadAtlasCKR <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                             LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    # load file
    pca.object <- read.table(PCA.OUTPUT, sep="\t", header=TRUE)
    pca.object <- pca.object %>% filter(Chain1 == "B" & Chain2 =="B")
    
    # substitute for BW and/or CCN names
    object.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.LIG.NAME][match(pca.object[ ,"ResNum1"], lookup.bwccn[ ,LOOK.LIG.RES])])
    colnames(object.comm.name) <- c("source_gnccn")
    ligand.comm.name <- as.data.frame(lookup.bwccn[ , LOOK.LIG.NAME][match(pca.object[ ,"ResNum2"], lookup.bwccn[ , LOOK.LIG.RES])])
    colnames(ligand.comm.name) <- c("target_gnccn")
    pca.object <- cbind(pca.object, object.comm.name)
    pca.object <- cbind(pca.object, ligand.comm.name)
    
    # add file name and PDB "class" (eg ck_ckr, cc_dimer, cxc_dimer)
    pdb.des <- as.data.frame(rep(PDBOBJ, nrow(pca.object)))
    colnames(pdb.des) <- c("file")
    pca.object <- cbind(pca.object, pdb.des)
    class.des <- as.data.frame(rep(CLASS, nrow(pca.object)))
    colnames(class.des) <- c("class")
    pca.object <- cbind(pca.object, class.des)
    
    # return final df, remove other objects
    return(pca.object)
    rm(pca.object, object.comm.name, ligand.comm.name, pca.object, pdb.des, class.des)
  }
  
  LoadAtlas <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                        LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
    
    inter <- LoadAtlasInter(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                            LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
    ck <- LoadAtlasCK(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                            LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
    ckr <- LoadAtlasCKR(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                            LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
    master <- rbind(inter, ck, ckr)
  }
  
##### 1: CALCULATE INTERMOLECULAR CONTACTS FROM CK-CKR COMPLEXES ###############
  
  # load CCN and BW conversion file
  lookup.bwccn <- read.csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv")
  
  # (1) CXCL8:CXCR1
  cna.1ilp <- LoadAtlas("01_structure_contacts/data/raw/1ILPCKCLEAN1538673510.txt",
                        "clean_1ilp_ck", "ccn_1ilp_ck",     # object resno, ccn
                        "clean_1ilp_ckr", "bw_1ilp_ckr",    # ligand resno, bw
                        "1ilp",                             # PDB ID
                        "soluble")                          # PDB "class"

  
  # (2) CCL11:CCR3
  cna.2mpm <- LoadAtlas("01_structure_contacts/data/raw/2MPMCKCLEAN1538684712.txt",
                          "clean_2mpm_ck", "ccn_2mpm_ck",   # object resno, ccn
                          "clean_2mpm_ckr", "bw_2mpm_ckr",  # ligand resno, bw
                          "2mpm",                           # PDB ID
                          "nmr")                            # PDB "class"
  
  # (3) CXCL12:CXCR4
  cna.2n55 <- LoadAtlas("01_structure_contacts/data/raw/2N55CKCLEAN1538684756.txt",
                          "clean_2n55_ck", "ccn_2n55_ck",     # object resno, ccn
                          "clean_2n55_ckr", "bw_2n55_ckr",  # ligand resno, bw
                          "2n55",                           # PDB ID
                          "soluble")                        # PDB "class"
  
  # (4) vMIPII:CXCR4
  cna.4rws <- LoadAtlas("01_structure_contacts/data/raw/4RWSCKCLEAN1538684823.txt",
                          "clean_4rws_ck", "ccn_4rws_ck",     # object resno, ccn
                          "clean_4rws_ckr", "bw_4rws_ckr",  # ligand resno, bw
                          "4rws",                           # PDB ID
                          "full")                           # PDB "class"
  
  # (5) CX3CL1:US28
  cna.4xt1 <- LoadAtlas("01_structure_contacts/data/raw/4XT1CKCLEAN1538684872.txt",
                          "clean_4xt1_ck", "ccn_4xt1_ck",     # object resno, ccn
                          "clean_4xt1_ckr", "bw_4xt1_ckr",  # ligand resno, bw
                          "4xt1",                           # PDB ID
                          "full")                           # PDB "class"
  
  # (6) CCL5[5P7]:CCR5
  cna.5uiw <- LoadAtlas("01_structure_contacts/data/raw/5UIWCKCLEAN1538684915.txt",
                          "clean_5uiw_ck", "ccn_5uiw_ck",     # object resno, ccn
                          "clean_5uiw_ckr", "bw_5uiw_ckr",  # ligand resno, bw
                          "5uiw",                           # PDB ID
                          "full")                           # PDB "class"

  # (7) CXCL12:CXCR4
  cna.2k05 <- LoadAtlas("01_structure_contacts/data/raw/2K05CLEAN1538684566.txt",
                          "clean_2k05_ck", "ccn_2k05_ck",   # object resno, ccn
                          "clean_2k05_ckr", "bw_2k05_ckr",  # ligand resno, bw
                          "2k05",                           # PDB ID
                          "soluble")                        # PDB "class"
  
  # (8) CCL5:CCR5
  cna.6fgp <- LoadAtlas("01_structure_contacts/data/raw/6FGPCLEAN1538688371.txt",
                        "clean_6fgp_ck", "ccn_6fgp_ck",     # object resno, ccn
                        "clean_6fgp_ckr", "bw_6fgp_ckr",    # ligand resno, bw
                        "6fgp",                             # PDB ID
                        "soluble")                          # PDB "class"
  
  # (9) CX3CL1.35:US28
  cna.5wb2 <- LoadAtlas("01_structure_contacts/data/raw/5WB2CLEAN1538684959.txt",
                        "clean_5wb2_ck", "ccn_5wb2_ck",     # object resno, ccn
                        "clean_5wb2_ckr", "bw_5wb2_ckr",    # ligand resno, bw
                        "5wb2",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (10) CCL20:CCR6
  cna.6wwz <- LoadAtlas("01_structure_contacts/data/raw/6WWZCLEAN1598531253.txt",
                        "clean_6wwz_ck", "ccn_6wwz_ck",     # object resno, ccn
                        "clean_6wwz_ckr", "bw_6wwz_ckr",    # ligand resno, bw
                        "6wwz",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (11) CXCL8:CXCR2
  cna.6lfo <- LoadAtlas("01_structure_contacts/data/raw/6LFOCLEAN1598542701.txt",
                        "clean_6lfo_ck", "ccn_6lfo_ck",     # object resno, ccn
                        "clean_6lfo_ckr", "bw_6lfo_ckr",    # ligand resno, bw
                        "6lfo",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (12) CXCL12:CXCR4 **MODEL**
  cna.ngo <- LoadAtlas("01_structure_contacts/data/raw/NGOMODELCLEAN1600451722.txt",
                        "clean_ngo_ck", "ccn_ngo_ck",     # object resno, ccn
                        "clean_ngo_ckr", "bw_ngo_ckr",    # ligand resno, bw
                        "ngo",                            # PDB ID
                        "full")                           # PDB "class"
  
  
  # bind all contacts into a single df (9 PDBs), then remove
  cna.master <- rbind(cna.1ilp, cna.6fgp, cna.2mpm, cna.2n55, cna.2k05, 
                      cna.5uiw, cna.4rws, cna.4xt1, cna.5wb2, cna.6lfo, cna.6wwz, cna.ngo)

  rm(cna.1ilp, cna.6fgp, cna.2mpm, cna.2n55, cna.2k05, 
                      cna.5uiw, cna.4rws, cna.4xt1, cna.5wb2, cna.6lfo, cna.6wwz, cna.ngo)
  
  # add domain designations
  lookup.domain <- read.csv("01_structure_contacts/data/processed/lookup_gnccn_to_domain.csv")
  cna.master$dom1 <- lookup.domain$dom[match(unlist(cna.master[ ,"source_gnccn"]), lookup.domain$bwccn)]
  cna.master$dom2 <- lookup.domain$dom[match(unlist(cna.master[ ,"target_gnccn"]), lookup.domain$bwccn)]  
  rm(lookup.bwccn, lookup.domain)
  
  # remove unecessary columns
  cna.master <- cna.master %>% select(-PDB, -SS1, -SS2)
  
  # write df
  write_csv(cna.master, "01_structure_contacts/output/RIN_atom.csv")
  
  # get side-chain contacts only
  cna.master.sc <- cna.master %>% filter(Chain.Types == "S-S") %>% 
    select(-Number.of.atomic.contacts, -Atoms, -Chain.Types, -Distance) %>%
    unique()
  
  # make unique such that each residue can only make one unique contact
  cna.master <- cna.master %>% 
    select(-Number.of.atomic.contacts, -Atoms, -Chain.Types, -Distance) %>%
    unique()

  # write df
  write_csv(cna.master, "01_structure_contacts/output/RIN_residue.csv")
  write_csv(cna.master.sc, "01_structure_contacts/output/RIN_residue_sc.csv")
  