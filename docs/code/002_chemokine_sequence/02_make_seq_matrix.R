# Imports sequence alignments and converts to data frame. Adds sequence labels
# (e.g. CC versus CXC) and writes output *.csv file which can subsequently
# be used for classification (e.g. CC versus CXC classification and testing). 
# Note that KMAD algorithm introduced lower case letters where substitutions 
# were made, however these are recognized as different symbols during 
# classification. As such all symbols are made upper case via Align2DataFrame
# function. Note that in section (1) paralogs are imported and will be used
# for CC/CXC training and testing. Additionally, human chemokine paralogs and
# viral chemokines are imported in sections (2) and (3), respectively and 
# converted to data frames for assignment of position-specific prediction
# probability scores for all sequences contained in these two alignment sets 
# based on models trained with CC and CXC ortholog sequences from (1) 
#
# Note on viral chemokine sequences:
# (1) Combined 30 viral sequences with 46 human paralogs and aligned using 
#     ClustalO online: viral_human_ckrs_clustal.fasta
# (2) Cross checked following alignment pairs by comparing to structural
#     alignmnet: vMIP2-CCL2
# (3) Pasted 30 viral sequences plus CCL2 sequence into new alignment with
#     46 human paralogs from master alignment; used CCL2 as bridge to adjust
#     viral sequences to master alignment; ***insertions in viral sequences
#     relative to master alignment were deleted***
# (4) Final alignment: SEQUENCES_VIRAL_CK.fasta 

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) MAKE CK SEQ MATRIX LABELED WITH "CC" AND "CXC" ---------------------------
# import
data <- Align2DataFrame("data/sequence/chemokine/alignments/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.fasta", 
                        "data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, " | ", remove = TRUE) %>% dplyr::select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = case_when(
  grepl("cx3c" , data$protein) ~ "NA",
  grepl("cc" , data$protein) ~ "cc",
  grepl("cxc" , data$protein) ~ "cxc"
))

# filter & reorder
data <- data %>% filter(class == "cc" | class == "cxc")
data <- data %>% dplyr::select(protein, seq, class, NTc.Cm70:CT.308)

# write output
# write_csv(data, "data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv")
# LAST WRITTEN 20230924
rm(data)

# (2) MAKE "TEST SET" CONTAINING CHEMOKINE PARALOGS ----------------------------
# import
data <- Align2DataFrame("data/sequence/chemokine/alignments/ALL_para.fasta", 
                        "data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, " | ", remove = TRUE) %>% dplyr::select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = case_when(
  grepl("cx3c" , data$protein) ~ "non",
  grepl("cc" , data$protein) ~ "cc",
  grepl("cxc" , data$protein) ~ "cxc",
  protein == "xcl1" ~ "non",
  protein == "xcl2" ~ "non"
  
))

# reorder
data <- data %>% dplyr::select(protein, seq, class, NTc.Cm70:CT.308)

# write output
# write_csv(data, "data/sequence/chemokine/alignment_csv/ALL_para_df.csv")
# LAST WRITTEN 20230924
rm(data)

# (3)  MAKE "TEST SET" CONTAINING VIRAL CHEMOKINES -----------------------------
# import
data <- Align2DataFrame("data/sequence/chemokine/alignments/SEQUENCES_VIRAL_CK.fasta", 
                        "data/lookup/SEQUENCES_CK_ALL_ORTHO_PARA_NTERM_UPDATED.txt")

# add protein name
pname <- data %>% separate(seq, "|", remove = TRUE) %>% dplyr::select(1) 
colnames(pname) <- c("protein")
data <- cbind(pname, data)
rm(pname)

# cc/cxc classifier
data <- data %>% mutate(class = "non")

# reorder
data <- data %>% dplyr:select(protein, seq, class, NTc.Cm70:CT.308)

# write
# write_csv(data, "data/sequence/chemokine/alignment_csv/ALL_virus_df.csv")
# LAST WRITTEN 20230924
rm(data)

