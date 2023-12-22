# Imports sequence alignments and converts to data frame. Adds sequence labels
# (e.g. CC versus CXC) and writes output *.csv file which can subsequently
# be used for classification (e.g. CC versus CXC classification and testing). 
# Note that in section (1) paralogs are imported and will be used
# for CC/CXC training and testing. Additionally, human class A GPCRs 
# (including chemokine receptor paralogs) and viral chemokine receptors
# are imported in sections (2) and (3), respectively and converted to data 
# frames for assignment of position-specific prediction
# probability scores for all sequences contained in these two alignment sets 
# based on models trained with CC and CXC ortholog sequences from (1) 
#
# Note on viral chemokine receptors:
# (1) Combined 11 viral seqs with 23 human receptor paralogs and aligned using
#     ClustalO: viral_human_ckrs_clustal.fasta
# (2) GPCRdb has US28 sequence - compared Clustal alignment to GPCRdb alignment
#     between US28 and CXCR4; alignments match in core regions
# (3) Pasted 11 viral sequences + CXCR4 into new alignment with 23 human paralogs
#     from the final master alignment; used CXCR4 as "bridge" to adjust viral
#     sequences to match master alignment; ***insertions in viral sequences
#     relative to the master alignment were deleted***; gaps in C-terminus after 
#     NPxxY were remoevd and resiudes left aligned to TM7; ECL2 adjusted as with
#     class A
# (4) Final alignment called SEQUENCES_VIRAL_NTERM_ECL2_UPDATED.fasta 

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) MAKE CKR SEQ MATRIX LABELED WITH "CC" AND "CXC" --------------------------
# import
data <- Align2DataFrame("data/sequence/gpcr/alignments/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.fasta", 
                        "data/lookup/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, " | ", remove = TRUE) %>% select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = case_when(
  grepl("ccrl2" , data$protein) ~ "NA",
  grepl("cc" , data$protein) ~ "cc",
  grepl("cxc" , data$protein) ~ "cxc"
))

# filter & reorder
data <- data %>% filter(class == "cc" | class == "cxc")
data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)

# write output
write_csv(data, "data/sequence/gpcr/alignment_csv/ALL_cc_cxc_ortho_df.csv")
rm(data)


# (2) MAKE NEW "TEST SET" CONTAINING CLASS A GPCRS -----------------------------
# import
data <- Align2DataFrame("data/sequence/gpcr/alignments/ALL_classa.fasta", 
                        "data/lookup/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, "_", remove = TRUE) %>% select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = case_when(
  grepl("ccrl2" , data$protein) ~ "non",
  grepl("cckar" , data$protein) ~ "non",
  grepl("cc" , data$protein) ~ "cc",
  grepl("cxc" , data$protein) ~ "cxc"
))

data <- data %>% mutate(class = case_when(
  is.na(class) ~ "non",
  !is.na(class) ~ class
))

# reorder
data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)

# write
write_csv(data, "data/sequence/gpcr/alignment_csv/ALL_classa_df.csv")


# (3) MAKE 2nd "TEST SET" CONTAINING VIRAL GPCRS -------------------------------

# import
data <- Align2DataFrame("03_ckr_seq/data/processed/ALL_virus.fasta", 
                        "03_ckr_seq/data/raw/SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, "|", remove = TRUE) %>% select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = "non")

# reorder
data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)

# write
write_csv(data, "03_ckr_seq/data/processed/ALL_virus_df.csv")


# (4) MAKE SET CONTAIING ALL PARALOGS MINUS ACKRs ------------------------------
# import
data <- Align2DataFrame("data/sequence/gpcr/alignments/ALL_non_ackr_para.fasta", 
                        "data/lookup//SEQUENCES_CKR_ALL_ORTHO_PARA_NTERM_ECL2_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, "|", remove = TRUE) %>% select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = "non")

# reorder
data <- data %>% select(protein, seq, class, gnNTr.Cm50:gnCT.57)

# write
write_csv(data, "data/sequence/gpcr/alignment_csv/ALL_non_ackr_para_df.csv")
