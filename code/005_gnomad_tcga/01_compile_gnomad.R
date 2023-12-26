# Imports protein-specific Gnomad files (downloaded from Gnomad), compiles 
# missense variant frequency counts for all proteins 
# into *.csv file. Note that the following mismatches between Gnomad reference 
# and alignment were identified (6 total because of duplicates):
# [1] "Mismatch in CCL15 at 24 Gnomad: ILE Aln: THR"
# [1] "Mismatch in CCL15 at 24 Gnomad: ILE Aln: THR" <- duplicate
# [1] "Mismatch in CCL23 at 123 Gnomad: VAL Aln: MET"
# [1] "Mismatch in CCL23 at 123 Gnomad: VAL Aln: MET" <- duplicate
# [1] "Mismatch in CXCL12 at 93 Gnomad: MET Aln: GLU"
# [1] "Mismatch in CXCL12 at 93 Gnomad: MET Aln: GLU" <- duplicate
# [1] "Mismatch in CXCL12 at 92 Gnomad: LYS Aln: ARG"
# [1] "Mismatch in CXCL13 at 94 Gnomad: ARG Aln: SER"
# [1] "Mismatch in CXCL13 at 107 Gnomad: LYS Aln: ILE"
# No other mismatches identified including for GPCRs. No adjustments made 
# to correct for these

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) IMPORT SEQS, MAKE VECTORS, PARSE GNOMAD - CHEMOKINE ----------------------
# import
aln.ck <- read.fasta("data/sequence/chemokine/alignments/ALL_para.fasta", as.string = T)
aln.ck <- sapply(aln.ck, as.character)

# get chemokine GENE NAMES
name.ck <- toupper(str_split(names(aln.ck), "_", simplify=T)[, 1]) # get chemokine names
names(aln.ck) <- name.ck # set chemokine names
ck.table <- data.frame(gene_symbol = name.ck, aln = aln.ck, stringsAsFactors=F) # create gene name / sequence table

# KEY STEP 1 - FUNCTION 1 - GNOMAD to TABLE
ck.variants <- adply(ck.table, 1, GetAlignedVariantsCK) 

# KEY STEP 2 - FUNCTION 2 - GNOMAD to ALN MATRIX (no. muts at each pos)
ck.variant.matrix <- MakeVariantMatrix(aln.ck, ck.variants)
ck.variant.matrix.freq <- MakeVariantMatrix(aln.ck, ck.variants, allele.count =T)
ck.variant.matrix.freq2 <- MakeVariantMatrixFreq(aln.ck, ck.variants, allele.count =T)

# removed used objects
rm(aln.ck, name.ck) #  objects


# (2) IMPORT SEQS, MAKE VECTORS, PARSE GNOMAD - RECEPTOR -----------------------
# import
aln.ckr <- read.fasta("data/sequence/gpcr/alignments/ALL_para.fasta", as.string = T)
aln.ckr <- sapply(aln.ckr, as.character)

# get receptor GENE NAMES
name.ckr <- toupper(str_split(names(aln.ckr), "_", simplify=T)[, 1])
names(aln.ckr) <- name.ckr
ckr.table <- data.frame(gene_symbol = name.ckr, aln = aln.ckr, stringsAsFactors=F)

# KEY STEP 1 - FUNCTION 1 - GNOMAD to TABLE
ckr.variants <- adply(ckr.table, 1, GetAlignedVariantsCKR)

# KEY STEP 2 - FUNCTION 2 - GNOMAD to ALN MATRIX (no. muts at each pos)
ckr.variant.matrix <- MakeVariantMatrix(aln.ckr, ckr.variants)
ckr.variant.matrix.freq <- MakeVariantMatrix(aln.ckr, ckr.variants, allele.count = T)
ckr.variant.matrix.freq2 <- MakeVariantMatrixFreq(aln.ckr, ckr.variants, allele.count = T)

# removed used objects
rm(aln.ckr, name.ckr) 

# (3) TIDY AND WRITE OUTPUT ----------------------------------------------------
# convert to data frame
ck.variant.matrix <- as.data.frame(ck.variant.matrix)
ck.variant.matrix.freq <- as.data.frame(ck.variant.matrix.freq)
ck.variant.matrix.freq2 <- as.data.frame(ck.variant.matrix.freq2)
ckr.variant.matrix <- as.data.frame(ckr.variant.matrix)
ckr.variant.matrix.freq <- as.data.frame(ckr.variant.matrix.freq)
ckr.variant.matrix.freq2 <- as.data.frame(ckr.variant.matrix.freq2)

# add names
ccn <- c(read.table("data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
colnames(ck.variant.matrix) <- c(ccn)
colnames(ck.variant.matrix.freq) <- c(ccn)
colnames(ck.variant.matrix.freq2) <- c(ccn)

gn <- c(read.table("data/lookup/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ",", colClasses = "character"))
colnames(ckr.variant.matrix) <- c(gn)
colnames(ckr.variant.matrix.freq) <- c(gn)
colnames(ckr.variant.matrix.freq2) <- c(gn)

# remove objects
rm(gn, ccn)

# add chemokine/receptor names
ck.variant.matrix <- cbind(rownames(ck.variant.matrix), ck.variant.matrix)
colnames(ck.variant.matrix)[1] <- c("protein")
ckr.variant.matrix <- cbind(rownames(ckr.variant.matrix), ckr.variant.matrix)
colnames(ckr.variant.matrix)[1] <- c("protein")
ck.variant.matrix.freq <- cbind(rownames(ck.variant.matrix.freq), ck.variant.matrix.freq)
colnames(ck.variant.matrix.freq)[1] <- c("protein")
ckr.variant.matrix.freq <- cbind(rownames(ckr.variant.matrix.freq), ckr.variant.matrix.freq)
colnames(ckr.variant.matrix.freq)[1] <- c("protein")
ck.variant.matrix.freq2 <- cbind(rownames(ck.variant.matrix.freq2), ck.variant.matrix.freq2)
colnames(ck.variant.matrix.freq2)[1] <- c("protein")
ckr.variant.matrix.freq2 <- cbind(rownames(ckr.variant.matrix.freq2), ckr.variant.matrix.freq2)
colnames(ckr.variant.matrix.freq2)[1] <- c("protein")

# gather
ck.variant.matrix <- ck.variant.matrix %>%
  gather(ccn, snp_count, 2:ncol(ck.variant.matrix))
ck.variant.matrix.freq <- ck.variant.matrix.freq %>%
  gather(ccn, snp_freq_count, 2:ncol(ck.variant.matrix.freq))
ck.variant.matrix.freq2 <- ck.variant.matrix.freq2 %>%
  gather(ccn, snp_freq_count, 2:ncol(ck.variant.matrix.freq2))
ckr.variant.matrix <- ckr.variant.matrix %>%
  gather(gn, snp_count, 2:ncol(ckr.variant.matrix))
ckr.variant.matrix.freq <- ckr.variant.matrix.freq %>%
  gather(gn, snp_freq_count, 2:ncol(ckr.variant.matrix.freq))
ckr.variant.matrix.freq2 <- ckr.variant.matrix.freq2 %>%
  gather(gn, snp_freq_count, 2:ncol(ckr.variant.matrix.freq2))

# join SNP counts and freq
ck.variant.matrix <- left_join(ck.variant.matrix, ck.variant.matrix.freq)
ckr.variant.matrix <- left_join(ckr.variant.matrix, ckr.variant.matrix.freq)
rm(ck.variant.matrix.freq, ckr.variant.matrix.freq)

# for table version add gn/ccn
lookup <- t(read.table("data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ","))
lookup <- as.data.frame(lookup)
colnames(lookup) <- c("ccn")
lookup$aln_pos <- 1:nrow(lookup)
ck.variants <- left_join(ck.variants, lookup)

lookup <- t(read.table("data/lookup/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ","))
lookup <- as.data.frame(lookup)
colnames(lookup) <- c("gn")
lookup$aln_pos <- 1:nrow(lookup)
ckr.variants <- left_join(ckr.variants, lookup)

# write output
# write_csv(ck.variant.matrix, "data/variant/gnomad/processed/CK_GNOMAD_COUNTS_FREQ.csv") # LAST WRITTEN 20231224
# write_csv(ck.variants, "data/variant/gnomad/processed/CK_GNOMAD_TABLE.csv") # LAST WRITTEN 20231224
# write_csv(ck.variant.matrix.freq2, "data/variant/gnomad/processed/CK_GNOMAD_COUNTS_FREQ_PCT.csv") # LAST WRITTEN 20231224

# write_csv(ckr.variant.matrix, "data/variant/gnomad/processed/CKR_GNOMAD_COUNTS_FREQ.csv") # LAST WRITTEN 20230927
# write_csv(ckr.variants, "data/variant/gnomad/processed/CKR_GNOMAD_TABLE.csv")
# write_csv(ckr.variant.matrix.freq2, "data/variant/gnomad/processed/CKR_GNOMAD_COUNTS_FREQ_PCT.csv") # LAST WRITTEN 20231204

# rm(ck.variants, ckr.variants, ck.variant.matrix, ckr.variant.matrix,
#    ck.table, ckr.table)
