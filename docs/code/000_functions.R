# Functions and descriptions listed below

Align2DataFrame <- function(ALIGNMENT, NAMES){
  
  # Description: imports alignment makes into data frame, adds common numbering
  # position names. Note that 
  #
  # Input:
  # ALIGNMENT = *.fasta file
  # NAMES = *.txt file containing character string of common numbers
  
  aln <- readAAMultipleAlignment(ALIGNMENT)
  aln.df <- as.matrix(aln)
  aln.df <- toupper(aln.df) # new 20202025
  seqname <- as.tibble(rownames(aln.df))
  colnames(seqname) <- c("seq")
  aln.df <- as.tibble(aln.df)
  gnccn <- t(read.table(NAMES, sep = ","))
  colnames(aln.df) <- gnccn
  aln.df <- cbind(seqname, aln.df)
  
  return(aln.df)
  rm(aln, aln.df, gnccn, seqname)
}

FormatVariantsGnomad <- function(var_table) {
  
  # define variants
  variant_types <- c("missense", "frameshift", "inframe insertion", "inframe deletion", "stop gained", "stop lost", "start lost")
  variant_colours <- c("blue", "indianred1",  "indianred2", "yellow", "red", "cyan2", "cyan3")
  names(variant_colours) <- variant_types
  
  var_subset <- 
    subset(var_table, !Annotation %in% c("5' UTR", "3' UTR", "synonymous", "intron",
                                         "non coding transcript exon", 
                                         "splice donor", "splice region",
                                         "splice acceptor", "upstream gene", 
                                         "mature miRNA", "downstream gene",
                                         "stop retained"))
  
  # filter low confidence
  var_subset$Flags[is.na(var_subset$Flags)] <- ""
  var_subset <- subset(var_subset, substr(Flags, 1, 2) != "LC")
  var_subset <- subset(var_subset, !(Annotation == "frameshift" & Protein.Consequence == ""))
  var_subset$MAF <- var_subset$Allele.Frequency
  var_subset$MAF[var_subset$MAF > 0.5] <- 1-var_subset$MAF[var_subset$MAF > 0.5]
  
  var_subset$MAF_cat <- factor(">=5%", levels=c("<1/10,000", "<1/1,000", "<1%", "<5%", ">=5%"))
  var_subset$MAF_cat[var_subset$MAF < 5/100] <- "<5%"
  var_subset$MAF_cat[var_subset$MAF < 1/100] <- "<1%"
  var_subset$MAF_cat[var_subset$MAF < 1/1000] <- "<1/1,000"
  var_subset$MAF_cat[var_subset$MAF < 1/10000] <- "<1/10,000"
  print("MAF categories assigned")
  
  # make the consequence names more pithy
  var_subset$Annotation[var_subset$Annotation == "missense_variant"] <- "missense"
  var_subset$Annotation[var_subset$Annotation == "frameshift_variant"] <- "frameshift"
  
  var_subset$Annotation <- factor(var_subset$Annotation, levels=variant_types)
  aa_consequence <- 
    str_match(var_subset$Protein.Consequence, "p.([A-Z][a-z][a-z])([0-9]+)(.+)")[, 2:4]
  aa_consequence <- 
    data.frame(aa_ref=aa_consequence[, 1], aa_pos=aa_consequence[, 2], aa_alt=aa_consequence[, 3], stringsAsFactors=F)
  var_subset[, c("aa_ref", "aa_pos", "aa_alt")] <- aa_consequence
  var_subset$aa_pos <- as.numeric(var_subset$aa_pos)
  
  var_subset
}

GetAlignedCancerVariants <- function(r) {
  gene_symbol <- r$gene_symbol
  print(gene_symbol)
  cancer_muts_subset <- subset(cancer.muts.parsed, gene_symbol == r$gene_symbol)
  
  if (!nrow(cancer_muts_subset)) {
    return(data.frame())
  }
  
  aln_split <- str_split(r$aln, "", simplify=T)
  gaps <- aln_split == "-"
  seq <- paste0(aln_split[!gaps])
  gap_offsets <- rep(NA, length(aln_split))
  gap_offsets[!gaps] <- 1:length(seq)
  
  print(subset(cancer_muts_subset, aa_pos > length(seq)))
  cancer_muts_subset <- subset(cancer_muts_subset, aa_pos <= length(seq))
  cbind(cancer_muts_subset, data.frame(gene_symbol=gene_symbol, aln_pos=sapply(cancer_muts_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
}

GetAlignedVariantsCK <- function(r) {
  gene_symbol <- r$gene_symbol
  print(gene_symbol)
  gnomad_file <- Sys.glob(paste0("data/variant/gnomad/raw/", gene_symbol, "_*.csv"))[1]
  print(gnomad_file)
  gnomad_vcf <- read.csv(gnomad_file, stringsAsFactors=F, na.strings="NA")
  gnomad_subset <- FormatVariantsGnomad(gnomad_vcf)
  # print(tail(gnomad_subset))
  
  aln_split <- str_split(r$aln, "", simplify=T)
  gaps <- aln_split == "-"
  seq <- paste0(aln_split[!gaps])
  gap_offsets <- rep(NA, length(aln_split))
  gap_offsets[!gaps] <- 1:length(seq)
  
  gnomad_subset <- subset(gnomad_subset, Annotation == "missense")
  print(subset(gnomad_subset, aa_pos > length(seq)))
  cbind(gnomad_subset, 
        data.frame(gene_symbol=gene_symbol, 
                   aln_pos=sapply(gnomad_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
}

GetAlignedVariantsCKR <- function(r) {
  gene_symbol <- r$gene_symbol
  print(gene_symbol)
  gnomad_file <- Sys.glob(paste0("data/variant/gnomad/raw/", gene_symbol, "_*.csv"))[1]
  print(gnomad_file)
  gnomad_vcf <- read.csv(gnomad_file, stringsAsFactors=F, na.strings="NA")
  gnomad_subset <- FormatVariantsGnomad(gnomad_vcf)
  # print(tail(gnomad_subset))
  
  aln_split <- str_split(r$aln, "", simplify=T)
  gaps <- aln_split == "-"
  seq <- paste0(aln_split[!gaps])
  gap_offsets <- rep(NA, length(aln_split))
  gap_offsets[!gaps] <- 1:length(seq)
  
  gnomad_subset <- subset(gnomad_subset, Annotation == "missense")
  print(subset(gnomad_subset, aa_pos > length(seq)))
  cbind(gnomad_subset, 
        data.frame(gene_symbol=gene_symbol, 
                   aln_pos=sapply(gnomad_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
}

GetContactList <- function(PDB1, PDB2){
  # Defines percentage of contacts shared (intersection) for 2 PDBs
  # compared to unique contacts in each PDB
  set1 <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
    filter(file == PDB1) %>%
    dplyr::select(file, source_gnccn, target_gnccn) %>%
    unique() %>%
    dplyr::select(-file)
  
  set2 <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
    filter(file == PDB2) %>%
    dplyr::select(file, source_gnccn, target_gnccn) %>%
    unique() %>%
    dplyr::select(-file)
  
  #--
  # set1 <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  #   filter(file == "8ic0") %>%
  #   dplyr::select(file, source_gnccn, target_gnccn) %>%
  #   unique() %>%
  #   dplyr::select(-file)
  # 
  # set2 <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  #   filter(file == "6lfo") %>%
  #   dplyr::select(file, source_gnccn, target_gnccn) %>%
  #   unique() %>%
  #   dplyr::select(-file)
  #--
  
  set1.not.set2 <- dplyr::setdiff(set1, set2)
  set1.not.set2$type <- c("a_not_b")
  set2.not.set1 <- dplyr::setdiff(set2, set1)
  set2.not.set1$type <- c("b_not_a")
  set1.i.set2 <- dplyr::intersect(set1, set2)
  set1.i.set2$type <- c("a_and_b")
  rm(set1, set2)
  df <- rbind(set1.not.set2, set2.not.set1, set1.i.set2)
  rm(set1.not.set2, set2.not.set1, set1.i.set2)
  return(df)
}

GetContactVenn <- function(PDB1, PDB2){
  # Defines percentage of contacts shared (intersection) for 2 PDBs
  # compared to unique contacts in each PDB
  set1 <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
    filter(file == PDB1) %>%
    dplyr::select(file, source_gnccn, target_gnccn) %>%
    unique() %>%
    dplyr::select(-file)
  
  set2 <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
    filter(file == PDB2) %>%
    dplyr::select(file, source_gnccn, target_gnccn) %>%
    unique() %>%
    dplyr::select(-file)
  
  set1.not.set2 <- dplyr::setdiff(set1, set2)
  set2.not.set1 <- dplyr::setdiff(set2, set1)
  set1.i.set2 <- dplyr::intersect(set1, set2)
  rm(set1, set2)
  
  return(nrow(set1.i.set2)/(nrow(set1.not.set2) + nrow(set2.not.set1) + nrow(set1.i.set2)))
}

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
  lookup.bwccn <- read.csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv")
  
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
  file1 <- strsplit(PDB1, split = '/')[[c(1,5)]]
  file1 <- strsplit(file1, split = '_')[[c(1,1)]]
  
  file2 <- strsplit(PDB2, split = '/')[[c(1,5)]]
  file2 <- strsplit(file2, split = '_')[[c(1,1)]]
  
  #--
  # file1 <- strsplit(pdbfiles[1], split = '/')[[c(1,4)]]
  # file1 <- strsplit(file1, split = '_')[[c(1,1)]]
  # 
  # file2 <- strsplit(pdbfiles[2], split = '/')[[c(1,4)]]
  # file2 <- strsplit(file2, split = '_')[[c(1,1)]]
  #--
  
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

GetPairwiseComplexCalphaRMSD <- function(PDB1, PDB2, 
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
  lookup.bwccn <- read.csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv")
  
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
  # select residues for structural alignment; will align based on receptor (chain B)
  # First need to identify shared positions among two PDB files using common numbering
  shared1 <- pdb1.atom %>% filter(chain == "B")
  shared1 <- shared1$gnccn
  shared2 <- pdb2.atom %>% filter(chain == "B")
  shared2 <- shared2$gnccn
  shared <- intersect(shared1, shared2)
  shared <- shared[!(shared %in% c("CT"))] # need to remove degenerate GNCCN designations
  rm(shared1, shared2)
  
  # now get resno and then indices for the set of shared positions
  shared1 <- pdb1.atom$resno[match(shared, pdb1.atom$gnccn)]
  shared2 <- pdb2.atom$resno[match(shared, pdb2.atom$gnccn)]
  
  pdb1.ind <- atom.select(pdb1, resno = c(shared1), chain = "B", elety='CA')
  pdb2.ind <- atom.select(pdb2, resno = c(shared2), chain = "B", elety='CA')
  
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
  # Now that the two PDBs are aligned based on the receptor position, you
  # want to calculate the RMSD of each pariwse CA of the chemokines. To do this
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
  file1 <- strsplit(PDB1, split = '/')[[c(1,5)]]
  file1 <- strsplit(file1, split = '_')[[c(1,1)]]
  
  file2 <- strsplit(PDB2, split = '/')[[c(1,5)]]
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

GetPairwiseReceptorCalphaRMSD <- function(PDB1, PDB2, 
                                          LOOK.CK1.GNCCN, LOOK.CK1.RESNO,
                                          LOOK.CKR1.GNCCN, LOOK.CKR1.RESNO,
                                          LOOK.CK2.GNCCN, LOOK.CK2.RESNO,
                                          LOOK.CKR2.GNCCN, LOOK.CKR2.RESNO){
  # Description:
  # Imports x2 PDB files, performs structural alignment of receptor, then 
  # calculates CA RMSD for all shared residues in receptor
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
  lookup.bwccn <- read.csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv")
  
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
  # will align based on receptor (chain B)
  # First need to identify shared positions among two PDB files using common numbering
  shared1 <- pdb1.atom %>% filter(chain == "B")
  shared1 <- shared1$gnccn
  shared2 <- pdb2.atom %>% filter(chain == "B")
  shared2 <- shared2$gnccn
  shared <- intersect(shared1, shared2)
  shared <- shared[!(shared %in% c("CT"))] # need to remove degenerate GNCCN designations
  rm(shared1, shared2)
  
  # now get resno and then indices for the set of shared positions
  shared1 <- pdb1.atom$resno[match(shared, pdb1.atom$gnccn)]
  shared2 <- pdb2.atom$resno[match(shared, pdb2.atom$gnccn)]
  
  pdb1.ind <- atom.select(pdb1, resno = c(shared1), chain = "B", elety='CA')
  pdb2.ind <- atom.select(pdb2, resno = c(shared2), chain = "B", elety='CA')
  
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
  # pariwse CA of the receptor To do this
  # you need to find common chemokine postions for which to compare CA RMSD
  
  # First need to identify shared positions among two PDB files using common numbering
  shared1 <- pdb1.atom %>% filter(chain == "B")
  shared1 <- shared1$gnccn
  shared2 <- pdb2.atom %>% filter(chain == "B")
  shared2 <- shared2$gnccn
  shared <- intersect(shared1, shared2)
  rm(shared1, shared2)
  
  # now get resno and then indices for the set of shared positions
  shared1 <- pdb1.atom$resno[match(shared, pdb1.atom$gnccn)]
  shared2 <- pdb2.atom$resno[match(shared, pdb2.atom$gnccn)]
  
  # create object to extract PDB ID from PDB path
  file1 <- strsplit(PDB1, split = '/')[[c(1,5)]]
  file1 <- strsplit(file1, split = '_')[[c(1,1)]]
  
  file2 <- strsplit(PDB2, split = '/')[[c(1,5)]]
  file2 <- strsplit(file2, split = '_')[[c(1,1)]]
  
  
  # make blank RMSD list
  rmsd.master <- as.data.frame(NULL)
  
  for (i in 1:length(shared)){
    pdb1.ind.i <- atom.select(pdb1, resno = c(shared1[i]),  chain = "B", elety='CA')
    pdb2.ind.i <- atom.select(pdb2, resno = c(shared2[i]),  chain = "B", elety='CA')
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

IndexToCodon <- function(idx) {
  # go from index to codon number and position within
  
  codon_n <- floor((idx-1) / 3) + 1
  codon_pos <- ((idx-1) %% 3) + 1
  
  return(c(codon_n, codon_pos))
}

GetPairwiseComplexCalphaRMSD <- function(PDB1, PDB2, 
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
  lookup.bwccn <- read.csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv")
  
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
  # select residues for structural alignment; will align based on receptor (chain B)
  # First need to identify shared positions among two PDB files using common numbering
  shared1 <- pdb1.atom %>% filter(chain == "B")
  shared1 <- shared1$gnccn
  shared2 <- pdb2.atom %>% filter(chain == "B")
  shared2 <- shared2$gnccn
  shared <- intersect(shared1, shared2)
  shared <- shared[!(shared %in% c("CT"))] # need to remove degenerate GNCCN designations
  rm(shared1, shared2)
  
  # now get resno and then indices for the set of shared positions
  shared1 <- pdb1.atom$resno[match(shared, pdb1.atom$gnccn)]
  shared2 <- pdb2.atom$resno[match(shared, pdb2.atom$gnccn)]
  
  pdb1.ind <- atom.select(pdb1, resno = c(shared1), chain = "B", elety='CA')
  pdb2.ind <- atom.select(pdb2, resno = c(shared2), chain = "B", elety='CA')
  
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
  # Now that the two PDBs are aligned based on the receptor position, you
  # want to calculate the RMSD of each pariwse CA of the chemokines. To do this
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
  file1 <- strsplit(PDB1, split = '/')[[c(1,5)]]
  file1 <- strsplit(file1, split = '_')[[c(1,1)]]
  
  file2 <- strsplit(PDB2, split = '/')[[c(1,5)]]
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

LoadAtlasInter <- function(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                           LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS){
  
  # Description: import Protein Contact Atlas (PCA) files and chemokine/GPCR common
  # numbering lookup tables, append common numbering to chemokine and GPCR
  # residue positions, add SSE information
  #
  # Input:
  # PCA.OUTPUT = .txt file from PCA
  # LOOK.OBJ.RES = column name from lookup table indicating chemokine residue  
  #     numbering for chemokine in question
  # LOOK.OBJ.NAME = column name from lookup table indicating CCN
  # LOOK.LIG.RES = column name from lookup table indicating GPCR residue  
  #     numbering for GPCR in question
  # LOOK.LIG.NAME = column name from lookup table indicating CGN
  # PDBOBJ = PDB ID
  # CLASS = "soluble" versus "full" complex
  
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
  
  # Description: See description associated with "LoadAtlasInter"
  
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
  
  # Description: See description associated with "LoadAtlasInter"
  

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
  
  # Description: See description associated with "LoadAtlasInter"
  
  inter <- LoadAtlasInter(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                          LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
  ck <- LoadAtlasCK(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                    LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
  ckr <- LoadAtlasCKR(PCA.OUTPUT, LOOK.OBJ.RES, LOOK.OBJ.NAME, 
                      LOOK.LIG.RES, LOOK.LIG.NAME, PDBOBJ, CLASS)
  master <- rbind(inter, ck, ckr)
}

MakeVariantMatrix <- function(aln, variants, allele.count=F, threeletter=T) {
  
  # (1) LOOP #1 - MAKE EMPTY MATRIX FOR MASTER CHEMOKINE ALIGNMENT
  # acts on "aln" which corresponds to ccl_seqs
  variant_matrix <- matrix(0, nrow=length(aln), ncol=nchar(aln[1])) # make blank matrix
  rownames(variant_matrix) <- names(aln) # change rownames
  
  for (s in rownames(variant_matrix)) {  # for each row...
    seq <- aln[s]                        # define new var that is a vector of 
    # seq positions from that row
    
    for (i in 1:(nchar(aln[1]))) {       # for each seq position...
      if (substr(seq, i, i) == '-') {    # replace dashes (ie no AA in alignment)
        variant_matrix[s, i] <- NA       # with NA
      }
    }
  }                                      # after this loop, all matrix positions
  # with zeros indicate that a chemokine
  # has a residue there, and all NAs
  # indicate that a chemokine has a dash
  
  # (2) LOOP #2 - MULTIPLE MANIPULATIONS OF EMPTY (ZERO) MATRIX
  # acts on "aln" which corresponds to ccl_seqs AND acts on ccl_variants which
  # is the master SNP table from GNOMAD
  
  for(i in 1:nrow(variants)) {          # for each row in GNOMAD table
    r <- variants[i, ]                  # make vector of values from table
    seq <- aln[r$gene_symbol]           # grab gene name (eg CCL1) and seq...
    
    
    # (2.1) IF-ELSE STATEMENT - IF "threeletter" IS TRUE
    # In both cases, fetches amino acid - outputs "K" for instance
    if (threeletter)
      aa_aln <- AMINO_ACID_CODE[toupper(substr(seq, r$aln_pos, r$aln_pos))] 
    else
      aa_aln <- toupper(substr(seq, r$aln_pos, r$aln_pos))
    
    
    # (2.2) IF STATEMENT - FIND MISMATCHES IN LISTED SNP AND REFERENCE
    if(toupper(r$aa_ref) != toupper(aa_aln)) {
      print(paste("Mismatch in", names(seq), "at", r$aa_pos, "Gnomad:", toupper(r$aa_ref), "Aln:", aa_aln))
    }
    
    
    # (2.3) IF-ELSE STATEMENT - IF "allele.count" IS TRUE
    if (allele.count) # if "allele.count" is TRUE...
      variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + r$Allele.Count
    # replace the matrix of 0's at a specific chemokine (row) AND 
    # alignment position (column) WITH existing value (ie 0) plus 
    # MAF (defined in function 1)
    else
      variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + 1
  }
  # otherwise add 1 (indicating frequency of 100% unless a SNP was listed)
  variant_matrix
}

MakeVariantMatrixFreq <- function(aln, variants, allele.count=F, threeletter=T) {
  
  # (1) LOOP #1 - MAKE EMPTY MATRIX FOR MASTER CHEMOKINE ALIGNMENT
  # acts on "aln" which corresponds to ccl_seqs
  variant_matrix <- matrix(0, nrow=length(aln), ncol=nchar(aln[1])) # make blank matrix
  rownames(variant_matrix) <- names(aln) # change rownames
  
  for (s in rownames(variant_matrix)) {  # for each row...
    seq <- aln[s]                        # define new var that is a vector of 
    # seq positions from that row
    
    for (i in 1:(nchar(aln[1]))) {       # for each seq position...
      if (substr(seq, i, i) == '-') {    # replace dashes (ie no AA in alignment)
        variant_matrix[s, i] <- NA       # with NA
      }
    }
  }                                      # after this loop, all matrix positions
  # with zeros indicate that a chemokine
  # has a residue there, and all NAs
  # indicate that a chemokine has a dash
  
  # (2) LOOP #2 - MULTIPLE MANIPULATIONS OF EMPTY (ZERO) MATRIX
  # acts on "aln" which corresponds to ccl_seqs AND acts on ccl_variants which
  # is the master SNP table from GNOMAD
  
  for(i in 1:nrow(variants)) {          # for each row in GNOMAD table
    r <- variants[i, ]                  # make vector of values from table
    seq <- aln[r$gene_symbol]           # grab gene name (eg CCL1) and seq...
    
    
    # (2.1) IF-ELSE STATEMENT - IF "threeletter" IS TRUE
    # In both cases, fetches amino acid - outputs "K" for instance
    if (threeletter)
      aa_aln <- AMINO_ACID_CODE[toupper(substr(seq, r$aln_pos, r$aln_pos))] 
    else
      aa_aln <- toupper(substr(seq, r$aln_pos, r$aln_pos))
    
    
    # (2.2) IF STATEMENT - FIND MISMATCHES IN LISTED SNP AND REFERENCE
    if(toupper(r$aa_ref) != toupper(aa_aln)) {
      print(paste("Mismatch in", names(seq), "at", r$aa_pos, "Gnomad:", toupper(r$aa_ref), "Aln:", aa_aln))
    }
    
    
    # (2.3) IF-ELSE STATEMENT - IF "allele.count" IS TRUE
    if (allele.count) # if "allele.count" is TRUE...
      variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + r$Allele.Frequency
    # replace the matrix of 0's at a specific chemokine (row) AND 
    # alignment position (column) WITH existing value (ie 0) plus 
    # MAF (defined in function 1)
    else
      variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + 1
  }
  # otherwise add 1 (indicating frequency of 100% unless a SNP was listed)
  variant_matrix
}

Raw2Motif <- function(FILE, LOOKUP_CLASS, MER, MASK){
  
  # import data
  data <- read_csv(FILE)
  data <- data[,-1]
  
  # extract chemokine name
  name <- strsplit(data$file, "/", fixed = TRUE)
  name <- t(as.data.frame(lapply(name, "[", 1)))
  colnames(name) <- ("protein")
  name <- strsplit(name[,1], "_", fixed = TRUE)
  name <- t(as.data.frame(lapply(name, "[", 1)))
  rownames(name) <- 1:nrow(name)
  colnames(name) <- c("protein")
  name <- as.data.frame(name)
  data <- cbind(data, name)
  rm(name)
  
  # add classification variables (CC, CXC, ACK, CX3L, XC)
  cc.cxc.ack <- read.csv(LOOKUP_CLASS)
  data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
  rm(cc.cxc.ack)
  
  # add mer designation
  data$mer <- c(MER)
  
  # reorder, filter, remove caps
  if(MER == "mer2"){
    
    data <- data %>% dplyr::select(protein, class, file, mer, motif_no, a, b) %>%
      filter(a != "-" & b != "-" )
    # force all caps
    a <- as.data.frame(toupper(data$a))
    b <- as.data.frame(toupper(data$b))
    colnames(a) <- c("a")
    colnames(b) <- c("b")
    data <- data %>% dplyr::select(-a,-b) %>% bind_cols(a,b) %>% 
      dplyr::select(protein, class, file, mer, motif_no, a, b)
    
  } else if (MER == "mer3"){
    
    data <- data %>% dplyr::select(protein, class, file, mer, motif_no, a, b, c) %>%
      filter(a != "-" & b != "-" & c != "-") 
    
    # force all caps
    a <- as.data.frame(toupper(data$a))
    b <- as.data.frame(toupper(data$b))
    c <- as.data.frame(toupper(data$c))
    colnames(a) <- c("a")
    colnames(b) <- c("b")
    colnames(c) <- c("c")
    data <- data %>% dplyr::select(-a,-b,-c) %>% bind_cols(a,b,c) %>% 
      dplyr::select(protein, class, file, mer, motif_no, a, b, c)
    
  } else if (MER == "mer4"){
    
    data <- data %>% dplyr::select(protein, class, file, mer, motif_no, a, b, c, d) %>%
      filter(a != "-" & b != "-" & c != "-" & d != "-") 
    
    # force all caps
    a <- as.data.frame(toupper(data$a))
    b <- as.data.frame(toupper(data$b))
    c <- as.data.frame(toupper(data$c))
    d <- as.data.frame(toupper(data$d))
    
    colnames(a) <- c("a")
    colnames(b) <- c("b")
    colnames(c) <- c("c")
    colnames(d) <- c("d")
    data <- data %>% dplyr::select(-a,-b,-c, -d) %>% bind_cols(a,b,c,d) %>% 
      dplyr::select(protein, class, file, mer, motif_no, a, b, c, d)
  }
  
  data$mask <- c(MASK)
  
  return(data)
  rm(data, a, b, c, d)
}

WriteCONECTcustom <- function(RINFILE, RINDF, PDBID, PDBFILE, OUTPUT){
  
  # import rinfile
  rin <- read_csv(RINFILE) %>% filter(file == PDBID)
  rin.df <- RINDF
  rin.df$sele <- c("yes")
  rin <- left_join(rin, rin.df)
  rin <- rin %>% filter(sele == "yes")
  rin <- rin %>% dplyr::select(-sele)
  
  # read PDB, make df, select relevant columns
  pdb <- read.pdb(PDBFILE)
  pdb_df <- as.data.frame(pdb$atom)
  pdb_conv <- pdb_df %>% dplyr::select(chain, resno, elety, eleno)
  pdb_conv <- pdb_conv %>% filter(elety == "CA")
  ck <- pdb_conv %>% filter(chain == "A")
  ckr <- pdb_conv %>% filter(chain == "B")
  
  # map atom indices to RIN file
  rin$ca1 <- ck$eleno[match(unlist(rin$ResNum1), ck$resno)]
  rin$ca2 <- ckr$eleno[match(unlist(rin$ResNum2), ckr$resno)]
  
  # clean up and write
  rin$CONECT <- c("CONECT")
  rin <- rin %>% dplyr::select(CONECT, ca1, ca2)
  write_csv(rin, OUTPUT)
  
  # return
  return(rin)
  
  # remove
  rm(rin, pdb, pdb_df, ck, ckr, pdb_conv)
}
