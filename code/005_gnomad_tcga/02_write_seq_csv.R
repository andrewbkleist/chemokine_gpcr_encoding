# Imports paralog alignments as *.fasta and writes output as 2-column *.csv
# with protein name (column 1) and aligned sequence (column 2) used for
# aligning TCGA mutation calls to master alignment in this paper

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CHEMOKINE ----------------------------------------------------------------
# import
aln.ck <- read.fasta("data/sequence/chemokine/alignments/ALL_para.fasta", as.string = T)
aln.ck <- sapply(aln.ck, as.character)

# get chemokine GENE NAMES
name.ck <- toupper(str_split(names(aln.ck), "_", simplify=T)[, 1]) # get chemokine names
names(aln.ck) <- name.ck # set chemokine names
ck.table <- data.frame(gene_symbol = name.ck, aln = aln.ck, stringsAsFactors=F) # create gene name / sequence table

# write output
# write_csv(ck.table, "data/sequence/chemokine/alignment_csv/CK_SEQ_TABLE.csv")
# WRITTEN 20231224

# (2) RECEPTOR -----------------------------------------------------------------
# import
aln.ckr <- read.fasta("data/sequence/gpcr/alignments/ALL_para.fasta", as.string = T)
aln.ckr <- sapply(aln.ckr, as.character)

# get receptor GENE NAMES
name.ckr <- toupper(str_split(names(aln.ckr), "_", simplify=T)[, 1])
names(aln.ckr) <- name.ckr
ckr.table <- data.frame(gene_symbol = name.ckr, aln = aln.ckr, stringsAsFactors=F)

# write output
# write_csv(ckr.table, "data/sequence/chemokine/alignment_csv/CKR_SEQ_TABLE.csv")
# WRITTEN 20230929

