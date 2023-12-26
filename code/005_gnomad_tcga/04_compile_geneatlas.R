source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

get_aligned_biobank <- function(r) {
  gene_symbol <- r$gene_symbol
  print(gene_symbol)
  biobank_subset <- subset(rbind(ccl_biobank, ccr_biobank), Gene == r$gene_symbol)
  if(nrow(biobank_subset) == 0)
    return (data.frame())
  
  aln_split <- str_split(r$aln, "", simplify=T)
  gaps <- aln_split == "-"
  seq <- paste0(aln_split[!gaps])
  gap_offsets <- rep(NA, length(aln_split))
  gap_offsets[!gaps] <- 1:length(seq)
  
  print(subset(biobank_subset, aa_pos > length(seq)))
  
  cbind(biobank_subset, 
        data.frame(gene_symbol=gene_symbol, 
                   aln_pos=sapply(biobank_subset$aa_pos, function(i) { which(i == gap_offsets) }), stringsAsFactors=F))
}

make_variant_matrix <- function(aln, variants, allele.count=F, threeletter=T) {
  
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


# Get a mapping of gene symbols to transcripts from Gnomad file names (elegant, I know)
id_map <- str_split(dir(path="data/variant/gnomad/raw/", pattern = "*.csv"), "_", simplify=T)
id_map <- as.data.frame(id_map)

ccl_biobank <- read.csv("data/variant/geneatlas/raw/chemokines_GA_hits_1e-2_processed_split copy.csv", sep=";", stringsAsFactors=F, header=T)
ccr_biobank <- read.csv("data/variant/geneatlas/raw/receptors_GA_hits_1e-2_processed_split copy.csv", sep=";", stringsAsFactors=F, header=T)

ccl_biobank <- ccl_biobank %>%
  mutate(Gene=toupper(Gene)) %>%
  filter(Transcript_ID %in% id_map$V2) %>%
  mutate(aa_ref=substr(AA_Consequences, 1, 1),
         aa_alt=substr(AA_Consequences, nchar(AA_Consequences), nchar(AA_Consequences)),
         aa_pos=as.numeric(substr(AA_Consequences, 2, nchar(AA_Consequences)-1)))

ccr_biobank <- ccr_biobank %>%
  mutate(Gene=toupper(Gene)) %>%
  filter(Transcript_ID %in% id_map$V2) %>%
  mutate(aa_ref=substr(AA_Consequences, 1, 1),
         aa_alt=substr(AA_Consequences, nchar(AA_Consequences), nchar(AA_Consequences)),
         aa_pos=as.numeric(substr(AA_Consequences, 2, nchar(AA_Consequences)-1)))

# new
ccl_table <- read_csv("data/variant/gnomad/processed/CK_GNOMAD_TABLE.csv") %>%
  dplyr::select(gene_symbol, aln)
ccr_table <- read_csv("data/variant/gnomad/processed/CKR_GNOMAD_TABLE.csv") %>%
  dplyr::select(gene_symbol, aln)

ccl_biobank_aln <- adply(ccl_table, 1, get_aligned_biobank)
ccr_biobank_aln <- adply(ccr_table, 1, get_aligned_biobank)

# map gn/ccn positions
lookup <- t(read.table("data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt", sep = ","))
lookup <- as.data.frame(lookup)
colnames(lookup) <- c("ccn")
lookup$aln_pos <- 1:nrow(lookup)
ccl_biobank_aln <- left_join(ccl_biobank_aln, lookup)

lookup <- t(read.table("data/lookup/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt", sep = ","))
lookup <- as.data.frame(lookup)
colnames(lookup) <- c("gn")
lookup$aln_pos <- 1:nrow(lookup)
ccr_biobank_aln <- left_join(ccr_biobank_aln, lookup)

# rm
rm(lookup,ccl_biobank, ccl_table, ccr_biobank, ccr_table, id_map )

# clean up
ccl_biobank_aln <- ccl_biobank_aln %>%
  dplyr::select(-aln) %>%
  unique()
ccr_biobank_aln <- ccr_biobank_aln %>%
  separate(col = gn, into = c("a", "gn"), sep = "gn", remove = FALSE) %>%
  dplyr::select(-a, -aln) %>%
  unique()
  

# write
# write_csv(ccl_biobank_aln, "data/variant/geneatlas/processed/CK_GENEATLAS.csv") # WRITTEN 2023.12.25
# write_csv(ccr_biobank_aln, "data/variant/geneatlas/processed/CKR_GENEATLAS.csv")

