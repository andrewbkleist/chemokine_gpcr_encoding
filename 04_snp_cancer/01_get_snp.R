# Name:     01_get_snp.R
# Updated:  20201030
# Author:   Andrew Kleist/Greg Slodkowicz

# packages, working directory
library(tidyverse)
library(plyr)
library(seqinr)
library(stringr)
library(Biostrings)
AMINO_ACID_CODE <- toupper(AMINO_ACID_CODE)
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")

# NOTE (20201030):
# TOPIC:DFs GENERATED
# ck.variants = lists variants (AA sub and position)
# ck.variant.matrix = for each chemokine (rows) by each positions (cols)
#                     gives number of variants for that chemokine/position
#                     in other words, how many UNIQUE MUTATIONS/SUBS
#                     happen in population irrespective of their frequency?
#                     (e.g. A->C, A->G, A->D = count of 3);
#                     ....coded in output as "snp_count"....
# ck.variant.matrix.freq = same as above *matrix except freq not counts;
#                          functionally you ADD the "allele.count" column
#                          from each of the different SNPs; e.g. if 
#                          A->C has 3 and A->G has 1 and A->D has 5 your
#                          "frequency" is 9; only mutations from the 
#                          reference sequence are counted; ie some listed SNPs
#                          might have F instead of A as their reference,
#                          so mutations of this F are not counted toward
#                          either the "snp_count" or "snp_count_freq"
# ...same for receptors...
#
#
# NOTE (20201101):
# TOPIC: MISMATCHED SEQS
# Note that the following mismatches between Gnomad reference and alignment were
# identified:
# [1] "Mismatch in CCL15 at 24 Gnomad: ILE Aln: THR"
# [1] "Mismatch in CCL15 at 24 Gnomad: ILE Aln: THR" <- duplicate
# [1] "Mismatch in CCL23 at 123 Gnomad: VAL Aln: MET"
# [1] "Mismatch in CCL23 at 123 Gnomad: VAL Aln: MET" <- duplicate
# [1] "Mismatch in CXCL12 at 93 Gnomad: MET Aln: GLU"
# [1] "Mismatch in CXCL12 at 93 Gnomad: MET Aln: GLU" <- duplicate
# [1] "Mismatch in CXCL12 at 92 Gnomad: LYS Aln: ARG"
# [1] "Mismatch in CXCL13 at 94 Gnomad: ARG Aln: SER"
# [1] "Mismatch in CXCL13 at 107 Gnomad: LYS Aln: ILE"
# No adjustments were made to correct for these; 6 total; none for receptors
#
#
# NOTE (20201101):
# TOPIC: SPOT CHECKS
# (1) To ensure that mappings are appropriately done with altered alignments from
# original analysis, the following steps were taken:
#   - For the following checked that sequential order of mutations per position 
#     was equivalent to that at the old alignment in the N-terminus;
#     also made spot checks in the core to verify that core mutations were
#     not changed from prior:
#   - Chemokines CCL4, CXCL5, CXCL17
#   - Receptors CCR4, CCR9, CXCR2
#
# (2) Also checked that both SNP FREQ and COUNTS correspond to the correct
# sums from the raw data dable from Gnomad for the following:
#   - Chemokines CCL1, CCL4
#   - Receptors CXCR1, CCR8

##### FUNCTIONS ################################################################
  
  # FUNCTION 1 -----------------------------------------------------------------
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
  
  # FUNCTION 2 -----------------------------------------------------------------
  GetAlignedVariants <- function(r) {
    gene_symbol <- r$gene_symbol
    print(gene_symbol)
    gnomad_file <- Sys.glob(paste0("04_snp_cancer/data/raw/", gene_symbol, "_*.csv"))[1]
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
  
  # FUNCTION 3 -----------------------------------------------------------------
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
      if(toupper(r$aa_ref) != aa_aln) {
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
  
  
##### 1: SNP TABLE AND MATRIX GENERATION #######################################

  # (1.1) IMPORT SEQS, MAKE VECTORS, PARSE GNOMAD - CHEMOKINE ------------------

    # import
    aln.ck <- read.fasta("02_ck_seq/data/processed/ALL_para.fasta", as.string = T)
    aln.ck <- sapply(aln.ck, as.character)
    
    # get chemokine GENE NAMES
    name.ck <- toupper(str_split(names(aln.ck), "_", simplify=T)[, 1]) # get chemokine names
    names(aln.ck) <- name.ck # set chemokine names
    ck.table <- data.frame(gene_symbol = name.ck, aln = aln.ck, stringsAsFactors=F) # create gene name / sequence table
    
    # KEY STEP 1 - FUNCTION 1 - GNOMAD to TABLE
    ck.variants <- adply(ck.table, 1, GetAlignedVariants) 

    # KEY STEP 2 - FUNCTION 2 - GNOMAD to ALN MATRIX (no. muts at each pos)
    ck.variant.matrix <- MakeVariantMatrix(aln.ck, ck.variants)
    ck.variant.matrix.freq <- MakeVariantMatrix(aln.ck, ck.variants, allele.count =T)
    
    # removed used objects
    rm(aln.ck, name.ck) #  objects
    
    
  # (1.2) IMPORT SEQS, MAKE VECTORS, PARSE GNOMAD - RECEPTOR -------------------
    
    # import
    aln.ckr <- read.fasta("03_ckr_seq/data/processed/ALL_para.fasta", as.string = T)
    aln.ckr <- sapply(aln.ckr, as.character)
    
    # get receptor GENE NAMES
    name.ckr <- toupper(str_split(names(aln.ckr), "_", simplify=T)[, 1])
    names(aln.ckr) <- name.ckr
    ckr.table <- data.frame(gene_symbol = name.ckr, aln = aln.ckr, stringsAsFactors=F)
    
    # KEY STEP 1 - FUNCTION 1 - GNOMAD to TABLE
    ckr.variants <- adply(ckr.table, 1, GetAlignedVariants)
    
    # KEY STEP 2 - FUNCTION 2 - GNOMAD to ALN MATRIX (no. muts at each pos)
    ckr.variant.matrix <- MakeVariantMatrix(aln.ckr, ckr.variants)
    ckr.variant.matrix.freq <- MakeVariantMatrix(aln.ckr, ckr.variants, allele.count = T)
    
    # removed used objects
    rm(aln.ckr, name.ckr) # chemokine objects
    
    rm(AMINO_ACID_CODE, FormatVariantsGnomad, GetAlignedVariants, MakeVariantMatrix)
    
    
  # (1.3) TIDY AND WRITE OUTPUT ------------------------------------------------
    
    # convert to data frame
    ck.variant.matrix <- as.data.frame(ck.variant.matrix)
    ck.variant.matrix.freq <- as.data.frame(ck.variant.matrix.freq)
    ckr.variant.matrix <- as.data.frame(ckr.variant.matrix)
    ckr.variant.matrix.freq <- as.data.frame(ckr.variant.matrix.freq)
    
    # add names
    ccn <- c(read.table("02_ck_seq/data/raw/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ",", colClasses = "character"))
    colnames(ck.variant.matrix) <- c(ccn)
    colnames(ck.variant.matrix.freq) <- c(ccn)
    gn <- c(read.table("03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ",", colClasses = "character"))
    colnames(ckr.variant.matrix) <- c(gn)
    colnames(ckr.variant.matrix.freq) <- c(gn)
    
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
    
    # gather
    ck.variant.matrix <- ck.variant.matrix %>%
      gather(ccn, snp_count, 2:ncol(ck.variant.matrix))
    ck.variant.matrix.freq <- ck.variant.matrix.freq %>%
      gather(ccn, snp_freq_count, 2:ncol(ck.variant.matrix.freq))
    ckr.variant.matrix <- ckr.variant.matrix %>%
      gather(gn, snp_count, 2:ncol(ckr.variant.matrix))
    ckr.variant.matrix.freq <- ckr.variant.matrix.freq %>%
      gather(gn, snp_freq_count, 2:ncol(ckr.variant.matrix.freq))
    
    # join SNP counts and freq
    ck.variant.matrix <- left_join(ck.variant.matrix, ck.variant.matrix.freq)
    ckr.variant.matrix <- left_join(ckr.variant.matrix, ckr.variant.matrix.freq)
    rm(ck.variant.matrix.freq, ckr.variant.matrix.freq)
    
    # write output
    write_csv(ck.table, "04_snp_cancer/data/processed/CK_SEQ_TABLE.csv")
    # write_csv(ck.variants, "04_snp_cancer/output/CK_GNOMAD_TABLE.csv")
    write_csv(ck.variant.matrix, "04_snp_cancer/output/CK_GNOMAD_COUNTS_FREQ.csv")

    write_csv(ckr.table, "04_snp_cancer/data/processed/CKR_SEQ_TABLE.csv")
    # write_csv(ckr.variants, "04_snp_cancer/output/CKR_GNOMAD_TABLE.csv")
    write_csv(ckr.variant.matrix, "04_snp_cancer/output/CKR_GNOMAD_COUNTS_FREQ.csv")

    rm(ck.variants, ckr.variants, ck.variant.matrix, ckr.variant.matrix,
       ck.table, ckr.table)

    