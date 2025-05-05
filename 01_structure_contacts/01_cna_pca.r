# Name:     01_cna_pca.r
# Updated:  20220731
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE 1 (20200918): Complex preparation
# all complexes were prepared by (1) removing hydrogens (via Bio3d clean.pdb);
# (2) the chemokine was made "Chain A" and the receptor "Chain B";
# (3) The "cleaned" PDB was submitted to Protein Contact Atlas 
# (https:/www.mrc-lmb.cam.ac.uk/rajini/index.html) with default settings
# GPCRdb nomenclature (GN) and Common Chemokine Numbering (CCN) designations
# were mapped using manually compiled "lookup table"; GN positions were modified
# to be inclusive of the conserved N-terminal Cys (designated 1x22);
# "Normal" GN positions from GPCRdb are identical past position 1x30

# NOTE 2 (20200918): Spot checks
# Contact mappings were validated on 20200918 by spot checks...
# Chose 3 pairwise contacts - made sure spatially made sense & GN/CCN match
#     - 1ilp / 2k05 / 2mpm / 2n55 / 6fgp
#     - 4rws / 5uiw / 5wb2 / 6lfo / 6wwz / ngo model
#     - additional spot checks 20220921...
#     - 7xbx / 7vl9 / 7sk3 / 7o7f / 7f1t
#     - 7f1r / 6meo / 7xa3
# The following source/target/dom1/dom2 fields have NAs: 
# (1) 1ilp NH2 and ACA "atoms"
# (2) "GS" at the beginning of CXCR4 from 2n55

# NOTE 3 (20220731): CCL5:CCR5 complexes
# There are multiple CCL5:CCR5 variant complexes including (1) CCL5[5P7]:CCR5 
# (5UIW),(2) CCL5[6P4]:CCR5 (7O7F), (3) CCL5:CCR5 (7F1R), (4) CCL5:CCR5 (Zheng);
# (the last beig a model). The native complex (7F1R) is missig a large part of
# CCR5 TM1/-term and CCL5 N-loop so "misses" a large chunk of the native "site 1/2"
# interaction. Nevertheless all (model included) show conserved CK orientation
# and similar CK N-term positions through NTc.Cm3 as well as CKR N-term course;
# Note also the native complex has a CK beta-1-to-ECL2 dislufide; 
# Since each complex is a different chemokine:GPCR complex (except native vs. model)
# will include all at this point. Model included b/c it recapitulates "invisible"
# parts of the solved complex which omits density for key regions.

# NOTE 4 (20220731): CX3CL1 and complexes
# CX3CL1:US28 (4XT1) and CX3CL1.35:US28 (5WB2) are both included since they 
# represent different N-terminal variants of CX3CL1 with the same receptor.
# CX3CL1:US28:G (7RKN/7RKF/7RKM) represent other versions of the same CK:CKR complex
# (CX3CL1:US28) so these were excluded from this analysis (4XT1 is higher
# resolution and was chosen for this reason). CX3CL1:CX3CR1 (7XBX) is its own
# unique pairing and was included. Note that in this complex the CKR N-term
# is missing and part of the CK B3/70s loop is missing, so it too suffers 
# from missing contacts (mostly "site 1").

# NOTE 5 (20220731): CCL15:CCR1 complexes
# There are 2 CCL15:CCR1 complexes (7VL9 and 7VLA) but the ligand is the same
# (CCL15) just a truncation variant of itself. As such, these are considered
# the same ligand and the longer N-term version (26-92) was chosen.

# NOTE 6 (20220731): CCL3:CCR5 complexes
# There are 2 CCL3:CCR5 complexes (7F1T and 7F1Q) but the ligand is the same
# (CCL15) just without or with G-protein. 7F1T was chosen as higher resolution.


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
  lookup.bwccn <- read.csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20220729.csv")
  
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
  cna.ngo <- LoadAtlas("01_structure_contacts/data/raw/NGOMODELCLEAN1614894949_1A.txt",
                       "clean_ngo_ck", "ccn_ngo_ck",     # object resno, ccn
                       "clean_ngo_ckr", "bw_ngo_ckr",    # ligand resno, bw
                       "ngo",                            # PDB ID
                       "full")                           # PDB "class"
  
  # (13) CCL5:CCR5 **MODEL**
  cna.zheng <- LoadAtlas("01_structure_contacts/data/raw/ZHENGMODELCLEAN1614894991_1A.txt",
                         "clean_zheng_ck", "ccn_zheng_ck",     # object resno, ccn
                         "clean_zheng_ckr", "bw_zheng_ckr",    # ligand resno, bw
                         "zheng",                              # PDB ID
                         "full")                              # PDB "class"
  
  # (14) CX3CL1:CX3CR1
  cna.7xbx <- LoadAtlas("01_structure_contacts/data/raw/7XBXCLEAN1659321025.txt",
                        "clean_7xbx_ck", "ccn_7xbx_ck",     # object resno, ccn
                        "clean_7xbx_ckr", "bw_7xbx_ckr",    # ligand resno, bw
                        "7xbx",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (15) CCL15:CCR1
  cna.7vl9 <- LoadAtlas("01_structure_contacts/data/raw/7VL9CLEAN1659321185.txt",
                        "clean_7vl9_ck", "ccn_7vl9_ck",     # object resno, ccn
                        "clean_7vl9_ckr", "bw_7vl9_ckr",    # ligand resno, bw
                        "7vl9",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (16) CXCL12:ACKR3
  cna.7sk3 <- LoadAtlas("01_structure_contacts/data/raw/7SK3CLEAN1659321603.txt",
                        "clean_7sk3_ck", "ccn_7sk3_ck",     # object resno, ccn
                        "clean_7sk3_ckr", "bw_7sk3_ckr",    # ligand resno, bw
                        "7sk3",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (17) CCL5[6P4]:CCR5
  cna.7o7f <- LoadAtlas("01_structure_contacts/data/raw/7O7FCLEAN1659321637.txt",
                        "clean_7o7f_ck", "ccn_7o7f_ck",     # object resno, ccn
                        "clean_7o7f_ckr", "bw_7o7f_ckr",    # ligand resno, bw
                        "7o7f",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (18) CCL3:CCR5
  cna.7f1t <- LoadAtlas("01_structure_contacts/data/raw/7F1TCLEAN1659321674.txt",
                        "clean_7f1t_ck", "ccn_7f1t_ck",     # object resno, ccn
                        "clean_7f1t_ckr", "bw_7f1t_ckr",    # ligand resno, bw
                        "7f1t",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (19) CCL5:CCR5
  cna.7f1r <- LoadAtlas("01_structure_contacts/data/raw/7F1RCLEAN1659321718.txt",
                        "clean_7f1r_ck", "ccn_7f1r_ck",     # object resno, ccn
                        "clean_7f1r_ckr", "bw_7f1r_ckr",    # ligand resno, bw
                        "7f1r",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (20) gp120:CCR5
  cna.6meo <- LoadAtlas("01_structure_contacts/data/raw/6MEOCLEAN1659334429.txt",
                        "clean_6meo_ck", "ccn_6meo_ck",     # object resno, ccn
                        "clean_6meo_ckr", "bw_6meo_ckr",    # ligand resno, bw
                        "6meo",                             # PDB ID
                        "full")                             # PDB "class"
  
  # (21) CCL2:CCR2
  cna.7xa3 <- LoadAtlas("01_structure_contacts/data/raw/7XA3CLEAN1663721528.txt",
                        "clean_7xa3_ck", "ccn_7xa3_ck",     # object resno, ccn
                        "clean_7xa3_ckr", "bw_7xa3_ckr",    # ligand resno, bw
                        "7xa3",                             # PDB ID
                        "full")                             # PDB "class"
  
  
  
  # bind all contacts into a single df (9 PDBs), then remove
  cna.master <- rbind(cna.1ilp, cna.2mpm, cna.2n55, cna.4rws, cna.4xt1, 
                      cna.5uiw, cna.2k05, cna.6fgp, cna.5wb2, cna.6wwz, 
                      cna.6lfo, cna.ngo, cna.zheng, cna.7xbx, cna.7vl9, 
                      cna.7sk3, cna.7o7f, cna.7f1t, cna.7f1r, cna.6meo, 
                      cna.7xa3
                      )

  rm(cna.1ilp, cna.2mpm, cna.2n55, cna.4rws, cna.4xt1, 
     cna.5uiw, cna.2k05, cna.6fgp, cna.5wb2, cna.6wwz, 
     cna.6lfo, cna.ngo, cna.zheng, cna.7xbx, cna.7vl9, 
     cna.7sk3, cna.7o7f, cna.7f1t, cna.7f1r, cna.6meo, 
     cna.7xa3)
  
  # add domain designations
  lookup.domain <- read.csv("01_structure_contacts/data/processed/lookup_gnccn_to_domain.csv")
  cna.master$dom1 <- lookup.domain$dom[match(unlist(cna.master[ ,"source_gnccn"]), lookup.domain$bwccn)]
  cna.master$dom2 <- lookup.domain$dom[match(unlist(cna.master[ ,"target_gnccn"]), lookup.domain$bwccn)]  
  rm(lookup.bwccn, lookup.domain)
  
  # remove unecessary columns
  cna.master <- cna.master %>% select(-PDB, -SS1, -SS2)
  
  # write df
  # write_csv(cna.master, "01_structure_contacts/output/RIN_atom.csv") # LAST WRITTEN 20221004
  # data <- read.csv("01_structure_contacts/output/RIN_atom.csv")
  
  # get side-chain contacts only
  # cna.master.sc <- cna.master %>% filter(Chain.Types == "S-S") %>% 
  #   select(-Number.of.atomic.contacts, -Atoms, -Chain.Types, -Distance) %>%
  #   unique()
  
  # make unique such that each residue can only make one unique contact
  cna.master <- cna.master %>%
    select(-Number.of.atomic.contacts, -Atoms, -Chain.Types, -Distance) %>%
    unique()

  # write df
  # write_csv(cna.master, "01_structure_contacts/output/RIN_residue.csv") # LAST WRITTEN 20221004
  # write_csv(cna.master.sc, "01_structure_contacts/output/RIN_residue_sc.csv")
  
  
  