# Name:     02_get_cancer.R
# Updated:  20201030
# Author:   Andrew Kleist

# packages, working directory
library(tidyverse)
library(plyr)
library(seqinr)
library(stringr)
library(Biostrings)
library(biomaRt)
AMINO_ACID_CODE <- toupper(AMINO_ACID_CODE)
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")

# NOTE (20201101):
# TOPIC: SPOT CHECKS
# (1) To ensure that mappings are appropriately done with altered alignments from
# original analysis, the following steps were taken:
#   - For CCL1 checked that sequential order of mutations per position was 
#     equivalent to that at the old alignment in the N-terminus;
#     also made spot checks in the core to verify that core mutations were
#     not changed from prior
#   - Repeated with chemokines CCL5, CXCL5, CXCL12
#   - Repeated with receptors CCR2, CCR5, ACKR3
#
# (2) Also checked that both MUTATION COUNTS correspond to the correct
# sums from the raw data dable from TCGA (via Danial McGrail) for the following:
#   - Chemokines CXCL1, CCL18
#   - Receptors ACKR4, CCR2

##### FUNCTIONS ################################################################
    
  # FUNCTION 1 -----------------------------------------------------------------
  IndexToCodon <- function(idx) {
    # go from index to codon number and position within
    
    codon_n <- floor((idx-1) / 3) + 1
    codon_pos <- ((idx-1) %% 3) + 1
    
    return(c(codon_n, codon_pos))
  }

  # FUNCTION 2 -----------------------------------------------------------------
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
  
  # FUNCTION 3 -----------------------------------------------------------------
  MakeVariantMatrix <- function(aln, variants, freq=F, threeletter=T) {
    variant_matrix <- matrix(0, nrow=length(aln), ncol=nchar(aln[1])) # make matrix of 0s
    rownames(variant_matrix) <- names(aln)
    
    # (1) LOOP #1 - MAKE EMPTY MATRIX FOR MASTER CHEMOKINE ALIGNMENT
    # acts on "aln" which corresponds to ccl_seqs
    for (s in rownames(variant_matrix)) {     # for each row...
      seq <- aln[s]                           # define new var that is a vector
      # seq positions from that row
      
      for (i in 1:(nchar(aln[1]))) {          # for each seq position
        if (substr(seq, i, i) == '-') {       # replace dashes (ie no AA in alignment)
          variant_matrix[s, i] <- NA          # with NA
        }
      }
    }                                         # after this loop, all matrix positions
    # with zeros indicate that a chemokine
    # has a residue there, and all NAs
    # indicate that a chemokine has a dash
    
    
    # (2) LOOP #2 - MULTIPLE MANIPULATIONS OF EMPTY (ZERO) MATRIX
    # acts on "aln" which corresponds to ccl_seqs AND acts on ccl_variants which
    # is the master SNP table from GNOMAD
    for(i in 1:nrow(variants)) {              # for each row in ccl_cancer_muts table
      r <- variants[i, ]                      # make vector of values from table
      seq <- aln[r$gene_symbol]               # grab gene name (eg CCL1) and seq
      
      # (2.1) IF-ELSE STATEMENT - IF "threeletter" IS TRUE
      # In both cases, fetches amino acid - outputs "K" for instance
      if (threeletter)
        aa_aln <- AMINO_ACID_CODE[toupper(substr(seq, r$aln_pos, r$aln_pos))]
      else
        aa_aln <- toupper(substr(seq, r$aln_pos, r$aln_pos))
      
      # (2.2) IF STATEMENT - FIND MISMATCHES IN LISTED SNP AND REFERENCE
      if(toupper(r$aa_ref) != aa_aln) {
        print(paste("Mismatch in", names(seq), "at", r$aa_pos, "Gnomad:", toupper(r$aa_ref), "Aln:", aa_aln))
      }
      
      # (2.3) IF-ELSE STATEMENT - IF "freq" IS TRUE
      if (freq) # if "freq" is TRUE...
        variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + r$MAF
      else
        variant_matrix[r$gene_symbol, r$aln_pos] <- variant_matrix[r$gene_symbol, r$aln_pos] + 1
    }
    
    variant_matrix
  }


##### 1: CANCER TABLE AND MATRIX GENERATION ####################################

  # (1.1) IMPORT FILES, SEQS FROM ENSEMBL --------------------------------------

    # import files, remove NAs, CXCR8
    cancer.muts <- read.csv("04_snp_cancer/data/raw/chemokinesMut5_20190207 copy 2.csv", stringsAsFactors = F)
    cancer.muts <- cancer.muts %>% filter(!is.na(Mut_new))
    cancer.muts <- cancer.muts %>% filter(Gene != "GPR35") # remove CXCR8
    
    # fetch sequences based on Ensembl IDs from Daniel's table
    ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="http://feb2014.archive.ensembl.org")
    coding.seqs <- getBM(attributes=c('hgnc_symbol', 'ensembl_gene_id', 'ensembl_peptide_id', "ensembl_transcript_id", "coding"),
                         filters='ensembl_transcript_id', values=unique(cancer.muts$ENST_new), mart=ensembl)
    rm(ensembl)
  
  # (1.2) SELECT & ANNOTATE NONSENSE MUTATIONS WITHIN DANIEL'S TABLE -----------
    
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


  # (1.3) TRANSLATE NUCLEOTIDE TO AMINO ACID AND ADD TO DANIEL'S TABLE ---------

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

    
  # (1.4) MAP ANDY'S ALIGNMENT TO PARSED SEQUENCE POSITIONS --------------------
    
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
    ck.table <- read_csv("04_snp_cancer/data/processed/CK_SEQ_TABLE.csv") # 
    ck.cancer.muts <- adply(ck.table, 1, GetAlignedCancerVariants) # match sequences
    
    # subset MUTATION TABLE for RECEPTORS, add ALIGNMENT POSITION
    ckr.table <- read_csv("04_snp_cancer/data/processed/CKR_SEQ_TABLE.csv")
    ckr.cancer.muts <- adply(ckr.table, 1, GetAlignedCancerVariants)
    
    # remove used objects
    rm(ck.table, ckr.table)

    # import CHEMOKINE sequences
    aln.ck <- read.fasta("02_ck_seq/data/processed/ALL_para.fasta", as.string = T)
    aln.ck <- sapply(aln.ck, as.character)
    name.ck <- toupper(str_split(names(aln.ck), "_", simplify=T)[, 1]) # get chemokine names
    names(aln.ck) <- name.ck
    
    # import RECEPTOR sequences
    aln.ckr <- read.fasta("03_ckr_seq/data/processed/ALL_para.fasta", as.string=T)
    aln.ckr <- sapply(aln.ckr, as.character)
    name.ckr <- toupper(str_split(names(aln.ckr), "_", simplify=T)[, 1])
    names(aln.ckr) <- name.ckr
    
  # (1.4) CONVERT TO MATRIX & ADD NAMES ----------------------------------------
    
    # CHEMOKINE - make, write matrix
    cancer.muts.ck.matrix <- MakeVariantMatrix(aln.ck, ck.cancer.muts, threeletter=F)
    ccn <- c(read.table("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
    colnames(cancer.muts.ck.matrix) <- c(ccn)
    cancer.muts.ck.matrix <- cbind(rownames(cancer.muts.ck.matrix), cancer.muts.ck.matrix)
    colnames(cancer.muts.ck.matrix)[1] <- c("protein")
    cancer.muts.ck.matrix <- as.data.frame(cancer.muts.ck.matrix)
    cancer.muts.ck.table <- cancer.muts.ck.matrix %>%
      gather(ccn, cancer_mut_count, 2:ncol(cancer.muts.ck.matrix))
    
    # RECEPTOR - make matrix
    cancer.muts.ckr.matrix <- MakeVariantMatrix(aln.ckr, ckr.cancer.muts, threeletter=F)
    gn <- c(read.table("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ",", colClasses = "character"))
    colnames(cancer.muts.ckr.matrix) <- c(gn)
    cancer.muts.ckr.matrix <- cbind(rownames(cancer.muts.ckr.matrix), cancer.muts.ckr.matrix)
    colnames(cancer.muts.ckr.matrix)[1] <- c("protein")
    cancer.muts.ckr.matrix <- as.data.frame(cancer.muts.ckr.matrix)
    cancer.muts.ckr.table <- cancer.muts.ckr.matrix %>%
      gather(gn, cancer_mut_count, 2:ncol(cancer.muts.ckr.matrix))
    
    # write
    write_csv(cancer.muts.ck.table, "04_snp_cancer/output/CK_TCGA_COUNTS.csv")
    write_csv(cancer.muts.ckr.table, "04_snp_cancer/output/CKR_TCGA_COUNTS.csv")
    
    # remove used objects
    rm(cancer.muts.ck.matrix, cancer.muts.ckr.matrix, cancer.muts.parsed,
       ck.cancer.muts, ccn, ckr.cancer.muts,
       gn, AMINO_ACID_CODE,
       GetAlignedCancerVariants, IndexToCodon, MakeVariantMatrix,
       aln.ck, aln.ckr, name.ck, name.ckr)
    