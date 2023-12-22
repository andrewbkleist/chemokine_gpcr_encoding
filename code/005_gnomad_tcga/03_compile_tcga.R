# Imports protein-specific TCGA files. Results in following mismatches:
# [1] "Mismatch in CXCL13 at 94 Gnomad: R Aln: S"
# [1] "Mismatch in CXCL13 at 102 Gnomad: P Aln: V"
# [1] "Mismatch in CXCL13 at 93 Gnomad: K Aln: R"
# [1] "Mismatch in CXCL13 at 106 Gnomad: R Aln: K"
# Mismactehs are chemokine-only, none found for receptors

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) IMPORT FILES, SEQS FROM ENSEMBL ------------------------------------------
# import files, remove NAs, GPR35 
cancer.muts <- read.csv("data/variant/tcga/raw/chemokinesMut5_20190207.csv", stringsAsFactors = F)
cancer.muts <- cancer.muts %>% filter(!is.na(Mut_new))
cancer.muts <- cancer.muts %>% filter(Gene != "GPR35") # low evidence as CKR

# fetch sequences based on Ensembl IDs from Daniel's table
ensembl <- useEnsembl(biomart="ensembl",
                      dataset="hsapiens_gene_ensembl",
                      host="https://grch37.ensembl.org")

coding.seqs <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 
                                  'ensembl_peptide_id', "ensembl_transcript_id", 
                                  "coding"),
                     filters='ensembl_transcript_id', 
                     values=unique(cancer.muts$ENST_new), 
                     mart=ensembl)
rm(ensembl)

# (2) SELECT & ANNOTATE NONSENSE MUTATIONS WITHIN DANIEL'S TABLE ---------------
# select missense mutations and isolate NUCLEOTIE, FROM, and TO (2600 mutations)
# Note that there are two syntaxes, so will treat individually and recombine
cancer.muts <- cancer.muts %>% filter(Type == "Missense_Mutation")
cancer.muts$temp_index <- c(1:nrow(cancer.muts))
muts.parsed <- cancer.muts %>% 
  filter(grepl("c.[ACGT][0-9]+[ACGT]", Mut_new) |
           grepl("c.[0-9]+[ACGT]>[ACGT]", Mut_new)) %>%
  dplyr::select(Mut_new, temp_index)
colnames(muts.parsed)[1] <- c("Mut")

# isolate two different syntaxes, get NUCLEOTIDE, FROM, TO, & recombine
synt.1 <- muts.parsed %>% filter(grepl("c.[ACGT][0-9]+[ACGT]", Mut)) # isolate
synt.2 <- muts.parsed %>% filter(grepl("c.[0-9]+[ACGT]>[ACGT]", Mut)) # isolate
synt.1.ann <- data.frame(str_match(synt.1$Mut, "c.([ACGT])([0-9]+)([ACGT])"), stringsAsFactors = F) # separate
synt.2.ann <- data.frame(str_match(synt.2$Mut, "c.([0-9]+)([ACGT])>([ACGT])"), stringsAsFactors = F) # separate
colnames(synt.1.ann) <- c("Mut_new","nt_from","pos", "nt_to") # annotate colnames
colnames(synt.2.ann) <- c("Mut_new", "pos", "nt_from", "nt_to") # annotate colnames
synt.1.ann$temp_index <- synt.1$temp_index # bring temporary index and NUC/FROM/TO together
synt.2.ann$temp_index <- synt.2$temp_index # bring temporary index and NUC/FROM/TO together
muts.parsed <- dplyr::bind_rows(synt.1.ann, synt.2.ann) # bring together fully annotated mutations from both syntaxes
cancer.muts <- left_join(cancer.muts, muts.parsed, by = "temp_index") # bring together into original table
cancer.muts$pos <- as.numeric(cancer.muts$pos)

# remove used variables
rm(muts.parsed, synt.1, synt.2, synt.1.ann, synt.2.ann)


# (3) TRANSLATE NUCLEOTIDE TO AMINO ACID AND ADD TO DANIEL'S TABLE -------------
# utilizes coding.seqs (SEQUENCE TABLE) and cancer_muts (mutation-annotated TCGA TABLE)
cancer.muts.parsed <- adply(cancer.muts, 1, function(r) {
  print(r$ENST_new)
  coding.seqs <- coding.seqs$coding[coding.seqs$ensembl_transcript_id == r$ENST_new]
  print(nchar(coding.seqs))
  coding_split <- substring(coding.seqs, seq(1, nchar(coding.seqs)-2, 3), seq(3, nchar(coding.seqs), 3))
  
  codon_coords <- IndexToCodon(r$pos)
  aa_pos <- codon_coords[1]
  codon_pos <- codon_coords[2]
  
  codon_from <- coding_split[aa_pos]
  codon_to <- codon_from
  substr(codon_to, codon_pos, codon_pos) <- r$nt_to
  
  print(paste(aa_pos, codon_from, codon_pos, r$nt_from))
  if(substr(codon_from, codon_pos, codon_pos) != r$nt_from) {
    print(paste("Mismatch in", r$gene_id, "at pos", aa_pos, codon_pos, "(expected", substr(codon_from, codon_pos, codon_pos), "saw", r$nt_from, ")"))
  }
  
  aa_ref <- GENETIC_CODE[codon_from]
  aa_to <- GENETIC_CODE[codon_to]
  
  data.frame(aa_pos=aa_pos, aa_ref=aa_ref, aa_to=aa_to, stringsAsFactors=F)
})

# rename column
cancer.muts.parsed$gene_symbol <- cancer.muts.parsed$Gene # required for get_aligned_cancer_variants

# remove used objects
rm(cancer.muts, coding.seqs)


# (4) MAP ANDY'S ALIGNMENT TO PARSED SEQUENCE POSITIONS ------------------------
# for dual-named chemokines/receptors, update table with conventional name
cancer.muts.parsed <- cancer.muts.parsed %>% mutate(Gene = case_when(
  Gene == "PF4" ~ "CXCL4",
  Gene == "PF4V1" ~ "CXCL4L1",
  Gene == "PPBP" ~ "CXCL7",
  Gene == "IL8" ~ "CXCL8",
  Gene == "CX3CR1" ~ "CX3C1",
  Gene == "DARC" ~ "ACKR1",
  Gene != "PF4" | Gene != "PF4V1" | Gene != "PPBP" | Gene != "IL8" | Gene != "CX3CR1" | Gene != "DARC" ~ Gene
))

# subset MUTATION TABLE for CHEMOKINES, add ALIGNMENT POSITION
ck.table <- read_csv("data/sequence/chemokine/alignment_csv/CK_SEQ_TABLE.csv") # 
ck.cancer.muts <- adply(ck.table, 1, GetAlignedCancerVariants) # match sequences

# subset MUTATION TABLE for RECEPTORS, add ALIGNMENT POSITION
ckr.table <- read_csv("data/sequence/chemokine/alignment_csv/CKR_SEQ_TABLE.csv")
ckr.cancer.muts <- adply(ckr.table, 1, GetAlignedCancerVariants)

# remove used objects
rm(ck.table, ckr.table)

# import CHEMOKINE sequences
aln.ck <- read.fasta("data/sequence/chemokine/alignments/ALL_para.fasta", as.string = T)
aln.ck <- sapply(aln.ck, as.character)
name.ck <- toupper(str_split(names(aln.ck), "_", simplify=T)[, 1]) # get chemokine names
names(aln.ck) <- name.ck

# import RECEPTOR sequences
aln.ckr <- read.fasta("data/sequence/gpcr/alignments/ALL_para.fasta", as.string=T)
aln.ckr <- sapply(aln.ckr, as.character)
name.ckr <- toupper(str_split(names(aln.ckr), "_", simplify=T)[, 1])
names(aln.ckr) <- name.ckr

# (5) CONVERT TO MATRIX & ADD NAMES --------------------------------------------
# CHEMOKINE - make, write matrix
cancer.muts.ck.matrix <- MakeVariantMatrix(aln.ck, ck.cancer.muts, threeletter=F)
ccn <- c(read.table("data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
colnames(cancer.muts.ck.matrix) <- c(ccn)
cancer.muts.ck.matrix <- cbind(rownames(cancer.muts.ck.matrix), cancer.muts.ck.matrix)
colnames(cancer.muts.ck.matrix)[1] <- c("protein")
cancer.muts.ck.matrix <- as.data.frame(cancer.muts.ck.matrix)
cancer.muts.ck.table <- cancer.muts.ck.matrix %>%
  gather(ccn, cancer_mut_count, 2:ncol(cancer.muts.ck.matrix))

# RECEPTOR - make matrix
cancer.muts.ckr.matrix <- MakeVariantMatrix(aln.ckr, ckr.cancer.muts, threeletter=F)
gn <- c(read.table("data/lookup/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ",", colClasses = "character"))
colnames(cancer.muts.ckr.matrix) <- c(gn)
cancer.muts.ckr.matrix <- cbind(rownames(cancer.muts.ckr.matrix), cancer.muts.ckr.matrix)
colnames(cancer.muts.ckr.matrix)[1] <- c("protein")
cancer.muts.ckr.matrix <- as.data.frame(cancer.muts.ckr.matrix)
cancer.muts.ckr.table <- cancer.muts.ckr.matrix %>%
  gather(gn, cancer_mut_count, 2:ncol(cancer.muts.ckr.matrix))

# write
# write_csv(cancer.muts.ck.table, "data/variant/tcga/processed/CK_TCGA_COUNTS.csv")
# write_csv(cancer.muts.ckr.table, "data/variant/tcga/processed/CKR_TCGA_COUNTS.csv")
# WRITTEN 20230929

# remove used objects
rm(cancer.muts.ck.matrix, cancer.muts.ckr.matrix, cancer.muts.parsed,
   ck.cancer.muts, ccn, ckr.cancer.muts,
   gn,
   GetAlignedCancerVariants, IndexToCodon, MakeVariantMatrix,
   aln.ck, aln.ckr, name.ck, name.ckr)

