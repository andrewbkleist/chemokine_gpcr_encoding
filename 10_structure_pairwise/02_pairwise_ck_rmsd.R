# Description:
# Calculates pairwise RMSD of all chemokines involved in complexes

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)

################################################################################

GetPairwiseChemokineCalphaRMSD <- function(PDB1, PDB2, 
                                  LOOK.CK1.GNCCN, LOOK.CK1.RESNO,
                                  LOOK.CKR1.GNCCN, LOOK.CKR1.RESNO,
                                  LOOK.CK2.GNCCN, LOOK.CK2.RESNO,
                                  LOOK.CKR2.GNCCN, LOOK.CKR2.RESNO){
  # Description:
  # Imports x2 PDB files, performs structural alignment of receptor, then 
  # calculates CA RMSD for all shared residues in chemokine.
  #
  # Args:
  #   PDB1 <- path to first PDB file
  #   PDB2 <- path to second PDB file (coordinates of this file are changed to match PDB1)
  #   LOOK.CK1.GNCCN <- column name containing CCN for first chemokine
  #   LOOK.CK1.RESNO <- column name containing residue numbers for first chemokine
  #   ... (remaining Args beginning with "LOOK" follow the same convention) ...
  #
  # Returns:
  #   Data frame with CCN numbering (column 1) CA RMSD (column 2), and paired
  #   PDBs (columns 3,4)
  
  # import pdbs
  pdb1 <- read.pdb(PDB1)
  #--
  # pdb1 <- read.pdb(pdbfiles[1])
  #--
  pdb1.atom <- pdb1$atom
  pdb1.atom <- pdb1.atom %>% filter(elety == "CA")
  pdb2 <- read.pdb(PDB2)
  #--
  # pdb2 <- read.pdb(pdbfiles[9])
  #--
  pdb2.atom <- pdb2$atom
  pdb2.atom <- pdb2.atom %>% filter(elety == "CA")
  
  # (2) COMMON NOMENCLATURE ----------------------------------------------------
  # add common nomenclatures for each residue
  # substitute for BW and/or CCN names
  lookup.bwccn <- read.csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20220729.csv")
  
  pdb1.A <- pdb1.atom %>% filter(chain == "A")
  pdb1.A$gnccn <- lookup.bwccn[ , LOOK.CK1.GNCCN][match(pdb1.A[ ,"resno"], lookup.bwccn[ ,LOOK.CK1.RESNO])]
  
  pdb1.B <- pdb1.atom %>% filter(chain == "B")
  pdb1.B$gnccn <- lookup.bwccn[ , LOOK.CKR1.GNCCN][match(pdb1.B[ ,"resno"], lookup.bwccn[ ,LOOK.CKR1.RESNO])]
  
  #--
  # pdb1.A <- pdb1.atom %>% filter(chain == "A")
  # pdb1.A$gnccn <- lookup.bwccn[ , lookup.ck.gnccn[1]][match(pdb1.A[ ,"resno"], lookup.bwccn[ ,lookup.ck.resno[1]])]
  # 
  # pdb1.B <- pdb1.atom %>% filter(chain == "B")
  # pdb1.B$gnccn <- lookup.bwccn[ , lookup.ckr.gnccn[1]][match(pdb1.B[ ,"resno"], lookup.bwccn[ ,lookup.ckr.resno[1]])]
  #--
  
  pdb1.atom <- rbind(pdb1.A, pdb1.B)
  rm(pdb1.A, pdb1.B)

  pdb2.A <- pdb2.atom %>% filter(chain == "A")
  pdb2.A$gnccn <- lookup.bwccn[ , LOOK.CK2.GNCCN][match(pdb2.A[ ,"resno"], lookup.bwccn[ , LOOK.CK2.RESNO])]
  
  pdb2.B <- pdb2.atom %>% filter(chain == "B")
  pdb2.B$gnccn <- lookup.bwccn[ , LOOK.CKR2.GNCCN][match(pdb2.B[ ,"resno"], lookup.bwccn[ , LOOK.CKR2.RESNO])]
  
  #--
  # pdb2.A <- pdb2.atom %>% filter(chain == "A")
  # pdb2.A$gnccn <- lookup.bwccn[ , lookup.ck.gnccn[9]][match(pdb2.A[ ,"resno"], lookup.bwccn[ ,lookup.ck.resno[9]])]
  # 
  # pdb2.B <- pdb2.atom %>% filter(chain == "B")
  # pdb2.B$gnccn <- lookup.bwccn[ , lookup.ckr.gnccn[9]][match(pdb2.B[ ,"resno"], lookup.bwccn[ ,lookup.ckr.resno[9]])]
  #--
  
  pdb2.atom <- rbind(pdb2.A, pdb2.B)
  rm(pdb2.A, pdb2.B)
  rm(lookup.bwccn)
  
  # (3) STRUCTURAL ALIGNMENT ---------------------------------------------------
  # select residues for structural alignment; 
  # will align based on chemokine (chain A)
  # First need to identify shared positions among two PDB files using common numbering
  shared1 <- pdb1.atom %>% filter(chain == "A")
  shared1 <- shared1$gnccn
  shared2 <- pdb2.atom %>% filter(chain == "A")
  shared2 <- shared2$gnccn
  shared <- intersect(shared1, shared2)
  shared <- shared[!(shared %in% c("CT"))] # need to remove degenerate GNCCN designations
  rm(shared1, shared2)
  
  # now get resno and then indices for the set of shared positions
  shared1 <- pdb1.atom$resno[match(shared, pdb1.atom$gnccn)]
  shared2 <- pdb2.atom$resno[match(shared, pdb2.atom$gnccn)]

  pdb1.ind <- atom.select(pdb1, resno = c(shared1), chain = "A", elety='CA')
  pdb2.ind <- atom.select(pdb2, resno = c(shared2), chain = "A", elety='CA')
  
  # now do CA superposition based on shared common positions but using residue numvers
  xyz <- fit.xyz(fixed = pdb1$xyz, 
                 mobile = pdb2$xyz,
                 fixed.inds = pdb1.ind$xyz, 
                 mobile.inds = pdb2.ind$xyz)
  # the xyz object is a matrix with 1 row (1 frame; would be multiple rows/frames
  # in case of MD simulation), and (example) 27,330 columns (xyz coords of all atoms from
  # the "mobile" or second listed PDB file). Note that these are the "new",
  # aligned xyz coords, so you can overwrite the xyz coords of the original
  # PDB object with these to make an aligned object
  
  # replace xyz coordinates of PDB2 with new, aligned XYZ coords
  pdb2$xyz <- xyz
  
  # sanity check - write output of moved second PDB file to show that 
  # alignment worked
  # write.pdb(pdb=pdb2, file = "moved.pdb")
  
  # (4) CALCULATE RMSD ---------------------------------------------------------
  # Now that the two PDBs are aligned , you want to calculate the RMSD of each 
  # pariwse CA of the chemokines. To do this
  # you need to find common chemokine postions for which to compare CA RMSD
  
  # First need to identify shared positions among two PDB files using common numbering
  shared1 <- pdb1.atom %>% filter(chain == "A")
  shared1 <- shared1$gnccn
  shared2 <- pdb2.atom %>% filter(chain == "A")
  shared2 <- shared2$gnccn
  shared <- intersect(shared1, shared2)
  rm(shared1, shared2)
  
  # now get resno and then indices for the set of shared positions
  shared1 <- pdb1.atom$resno[match(shared, pdb1.atom$gnccn)]
  shared2 <- pdb2.atom$resno[match(shared, pdb2.atom$gnccn)]

  # create object to extract PDB ID from PDB path
  file1 <- strsplit(PDB1, split = '/')[[c(1,4)]]
  file1 <- strsplit(file1, split = '_')[[c(1,1)]]
  
  file2 <- strsplit(PDB2, split = '/')[[c(1,4)]]
  file2 <- strsplit(file2, split = '_')[[c(1,1)]]
  
  
  # make blank RMSD list
  rmsd.master <- as.data.frame(NULL)
  
  for (i in 1:length(shared)){
    pdb1.ind.i <- atom.select(pdb1, resno = c(shared1[i]),  chain = "A", elety='CA')
    pdb2.ind.i <- atom.select(pdb2, resno = c(shared2[i]),  chain = "A", elety='CA')
    rmsd.res <- 
      data.frame(paste0(shared[i]), 
                 rmsd(a=pdb1$xyz[,pdb1.ind.i$xyz], 
                      b=pdb2$xyz[,pdb2.ind.i$xyz]),
                 file1, 
                 file2
      )
    colnames(rmsd.res) <- c("gnccn", "RMSD", "file1", "file2")
    # rmsd.res$resno <- as.numeric(rmsd.res$resno)
    rmsd.master <- rbind(rmsd.master, rmsd.res)
  }
  
  return(rmsd.master)
  
}
  
  


################################################################################

# define PDBs and lookup column names
pdbfiles <- c("01_structure_contacts/data/pdbs/5uiw_ck_clean.pdb",
              "01_structure_contacts/data/pdbs/7o7f_clean.pdb",
              "01_structure_contacts/data/pdbs/7f1r_clean.pdb",
              "01_structure_contacts/data/pdbs/zheng_model_clean.pdb",
              "01_structure_contacts/data/pdbs/7vl9_clean.pdb",
              "01_structure_contacts/data/pdbs/7xa3_clean.pdb",
              "01_structure_contacts/data/pdbs/7f1t_clean.pdb",
              "01_structure_contacts/data/pdbs/6wwz_clean.pdb",
              "01_structure_contacts/data/pdbs/6lfo_clean.pdb",
              "01_structure_contacts/data/pdbs/ngo_model_clean.pdb",
              "01_structure_contacts/data/pdbs/7sk3_clean.pdb",
              "01_structure_contacts/data/pdbs/7xbx_clean.pdb",
              "01_structure_contacts/data/pdbs/4rws_ck_clean.pdb",
              "01_structure_contacts/data/pdbs/4xt1_ck_clean.pdb",
              "01_structure_contacts/data/pdbs/5wb2_clean.pdb")
lookup.ck.gnccn <- c("ccn_5uiw_ck",
                     "ccn_7o7f_ck",
                     "ccn_7f1r_ck",
                     "ccn_zheng_ck",
                     "ccn_7vl9_ck",
                     "ccn_7xa3_ck",
                     "ccn_7f1t_ck",
                     "ccn_6wwz_ck",
                     "ccn_6lfo_ck",
                     "ccn_ngo_ck",
                     "ccn_7sk3_ck",
                     "ccn_7xbx_ck",
                     "ccn_4rws_ck",
                     "ccn_4xt1_ck",
                     "ccn_5wb2_ck")
lookup.ck.resno <- c("clean_5uiw_ck",
                     "clean_7o7f_ck",
                     "clean_7f1r_ck",
                     "clean_zheng_ck",
                     "clean_7vl9_ck",
                     "clean_7xa3_ck",
                     "clean_7f1t_ck",
                     "clean_6wwz_ck",
                     "clean_6lfo_ck",
                     "clean_ngo_ck",
                     "clean_7sk3_ck",
                     "clean_7xbx_ck",
                     "clean_4rws_ck",
                     "clean_4xt1_ck",
                     "clean_5wb2_ck")
lookup.ckr.gnccn <- c("bw_5uiw_ckr",
                      "bw_7o7f_ckr",
                      "bw_7f1r_ckr",
                      "bw_zheng_ckr",
                      "bw_7vl9_ckr",
                      "bw_7xa3_ckr",
                      "bw_7f1t_ckr",
                      "bw_6wwz_ckr",
                      "bw_6lfo_ckr",
                      "bw_ngo_ckr",
                      "bw_7sk3_ckr",
                      "bw_7xbx_ckr",
                      "bw_4rws_ckr",
                      "bw_4xt1_ckr",
                      "bw_5wb2_ckr")
lookup.ckr.resno <- c("clean_5uiw_ckr",
                      "clean_7o7f_ckr",
                      "clean_7f1r_ckr",
                      "clean_zheng_ckr",
                      "clean_7vl9_ckr",
                      "clean_7xa3_ckr",
                      "clean_7f1t_ckr",
                      "clean_6wwz_ckr",
                      "clean_6lfo_ckr",
                      "clean_ngo_ckr",
                      "clean_7sk3_ckr",
                      "clean_7xbx_ckr",
                      "clean_4rws_ckr",
                      "clean_4xt1_ckr",
                      "clean_5wb2_ckr")


# loop all-by-all comparisons of CA RMSD
data <- as.data.frame(NULL)
for (i in 1:length(pdbfiles)){
  for(j in 1:length(pdbfiles)){
    temp <- GetPairwiseChemokineCalphaRMSD(pdbfiles[i], pdbfiles[j],
                                  lookup.ck.gnccn[i], lookup.ck.resno[i],
                                  lookup.ckr.gnccn[i], lookup.ckr.resno[i],
                                  lookup.ck.gnccn[j], lookup.ck.resno[j],
                                  lookup.ckr.gnccn[j], lookup.ckr.resno[j])
    data <- rbind(data, temp)
  }
}
rm(temp)

# remove self-by-self comparisons
data <- data %>% filter(file1 != file2)

# remove duplicates
cols = c(1:4)
newdf = data[,cols]
for (i in 1:nrow(data)){
  newdf[i, ] = sort(data[i,cols])
}
data <- data[!duplicated(newdf),]

# write csv
write_csv(data, "10_structure_pairwise/output/pairwise_chemokine_rmsd.csv")

