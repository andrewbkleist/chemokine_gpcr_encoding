# Imports residue interaction network (RIN) data and annotates with (i) 
# conservation, (ii) subfamily scores, (iii) SNP data, and (iv) cancer data for
# each chemokine and GPCR residue position involved. For modified and non-human
# positions does not annotate. Adds structure statistics for each contact. See
# Methods for discussion of structure curation. Only annotates RINs for "full" 
# (i.e. non-soluble/NMR) complexes. The gp120-CCR5 complex was excluded from
# structure contact calculations.

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) Contact analysis  --------------------------------------------------------
# import rin
rin <- read_csv("data/structure/processed/RIN_residue.csv") %>%
  filter(class == "full") %>% filter(Chain1 != Chain2) %>%
  dplyr::select(file, source_gnccn, target_gnccn, dom1, dom2)

# add representation of  each contact among *all* complexes
n.pdb <- rin %>% 
  filter(!(file %in% c("6meo"))) %>% # remove gp120:CCR5 complex
  dplyr::select(source_gnccn, target_gnccn) %>% 
  dplyr::count(source_gnccn, target_gnccn) 
colnames(n.pdb)[3] <- c("no_pdb")
rin <- left_join(rin, n.pdb)
rm(n.pdb)

# add number of unique residues contacted by each chemokine/GPCR residue WITHIN COMPLEX
n.resi.ck <- rin %>% dplyr::select(file, source_gnccn) %>% dplyr::count(file, source_gnccn)
colnames(n.resi.ck)[3] <- c("no_contacts_file_ck")
n.resi.ckr <- rin %>% dplyr::select(file, target_gnccn) %>% dplyr::count(file, target_gnccn)
colnames(n.resi.ckr)[3] <- c("no_contacts_file_ckr")
rin <- left_join(rin, n.resi.ck)
rin <- left_join(rin, n.resi.ckr)
rm(n.resi.ck, n.resi.ckr)

# add number of unique residues contacted by each chemokine/GPCR residue ACROSS COMPLEXES
n.resi.ck <- rin %>% 
  filter(file != "6meo") %>%
  # remove gp120-CCR5 complex so as not to skew contacts
  dplyr::select(source_gnccn, target_gnccn) %>% unique() %>% 
  dplyr::select(source_gnccn) %>% dplyr::count(source_gnccn)
colnames(n.resi.ck)[2] <- c("no_contacts_all_ck")

n.resi.ckr <- rin %>% 
  filter(file != "6meo") %>%
  # remove gp120-CCR5 complex 
  dplyr::select(source_gnccn, target_gnccn) %>% unique() %>% 
  dplyr::select(target_gnccn) %>% dplyr::count(target_gnccn) 
colnames(n.resi.ckr)[2] <- c("no_contacts_all_ckr")

rin <- left_join(rin, n.resi.ck)
rin <- left_join(rin, n.resi.ckr)

# add common names so you can map conservation information
rin <- rin %>% mutate(ck = case_when(
  file == "5uiw" ~ "ccl5",
  file == "4rws" ~ "vmipii",
  file == "4xt1" ~ "cx3cl1",
  file == "5wb2" ~ "cx3cl1",
  file == "6lfo" ~ "cxcl8",
  file == "6wwz" ~ "ccl20",
  file == "ngo" ~ "cxcl12",
  file == "zheng" ~ "ccl5",
  file == "7xbx" ~ "cx3cl1",
  file == "7vl9" ~ "ccl15",
  file == "7sk3" ~ "cxcl12",
  file == "7o7f" ~ "ccl5",
  file == "7f1t" ~ "ccl3",
  file == "7f1r" ~ "ccl5",
  file == "6meo" ~ "gp120",
  file == "7xa3" ~ "ccl2",
  file == "8ic0" ~ "cxcl8"
)) %>% mutate(ckr = case_when(
  file == "5uiw" ~ "ccr5",
  file == "4rws" ~ "cxcr4",
  file == "4xt1" ~ "us28",
  file == "5wb2" ~ "us28",
  file == "6lfo" ~ "cxcr2",
  file == "6wwz" ~ "ccr6",
  file == "ngo" ~ "cxcr4",
  file == "zheng" ~ "ccr5",
  file == "7xbx" ~ "cx3cr1",
  file == "7vl9" ~ "ccr1",
  file == "7sk3" ~ "ackr3",
  file == "7o7f" ~ "ccr5",
  file == "7f1t" ~ "ccr5",
  file == "7f1r" ~ "ccr5",
  file == "6meo" ~ "ccr5",
  file == "7xa3" ~ "ccr2",
  file == "8ic0" ~ "cxcr1"
))
rin <- rin %>% dplyr::select(file, ck, ckr, source_gnccn, target_gnccn, dom1, dom2, 
                      no_pdb, no_contacts_file_ck, no_contacts_file_ckr,
                      no_contacts_all_ck, no_contacts_all_ckr)

rm(n.resi.ck, n.resi.ckr)


# (2) ADD CONSERVATION ---------------------------------------------------------
# import conservation chemokine
ck.cons  <- read_csv("data/integrated/CK_CONS_CCCXC_SNP_CAN.csv")
colnames(ck.cons) <- c("ck", "source_gnccn","dom1", "all_para_ck", "all_cc_cxc_para_ck",
                       "cc_para_ck", "cxc_para_ck", "ortho_cons_ck",  "ngaps_ck",  "pctgaps_ck",
                       "cc_cxc_lr_ck", "cc_cxc_lr_sd_ck", "cc_cxc_lr_score_ck", "cc_cxc_lr_score_sd_ck",
                       "snp_count_ck", "snp_freq_count_ck", "cancer_mut_count_ck")
ck.ortho.cons <- ck.cons %>% dplyr::select(ck, source_gnccn, ortho_cons_ck) %>% 
  separate(ck, sep = "_", into = "ck")
ck.cons <- ck.cons %>% 
  dplyr::select("source_gnccn","dom1", "all_para_ck",  "all_cc_cxc_para_ck", 
         "cc_para_ck", "cxc_para_ck", "cc_cxc_lr_ck", "cc_cxc_lr_sd_ck") %>%
  unique()

# import conservation receptor
ckr.cons  <- read_csv("data/integrated/CKR_CONS_CCCXC_SNP_CAN.csv")  
colnames(ckr.cons) <- c("ckr", "target_gnccn","dom2", "all_para_ckr", "all_non_ackr_para_ckr", "all_cc_cxc_para_ckr", "all_classa_ckr",  
                        "cc_para_ckr", "cxc_para_ckr", "ack_para_ckr", "ortho_cons_ckr",  "ngaps_ckr",  "pctgaps_ckr",
                        "cc_cxc_lr_ckr", "cc_cxc_lr_sd_ckr", "cc_cxc_lr_score_ckr", 
                        "cc_cxc_lr_score_sd_ckr",
                        "snp_count_ckr", "snp_freq_count_ckr", 
                        "cancer_mut_count_ckr", "resid")
ckr.cons <- ckr.cons %>% separate(col = target_gnccn,  into  = c("temp","target_gnccn"), sep = "gn", remove = FALSE)
ckr.cons <- ckr.cons  %>% dplyr::select(-temp)

ckr.ortho.cons <- ckr.cons %>% dplyr::select(ckr, target_gnccn, ortho_cons_ckr) %>% 
  separate(ckr, sep = "_", into = "ckr")
ckr.cons <- ckr.cons %>% 
  dplyr::select("target_gnccn","dom2", "all_para_ckr", "all_non_ackr_para_ckr", "all_cc_cxc_para_ckr","all_classa_ckr",  
         "cc_para_ckr", "cxc_para_ckr", "ack_para_ckr", "cc_cxc_lr_ckr", "cc_cxc_lr_sd_ckr") %>%
  unique()

#  combine generalized scores (ie protein-independent)
rin <- left_join(rin, ck.cons)
rin <- left_join(rin, ckr.cons)
ck.ortho.cons$ck <- tolower(ck.ortho.cons$ck)
ckr.ortho.cons$ckr <- tolower(ckr.ortho.cons$ckr)


# combine protein-specific ortholog scores CHEMOKINE
rin <- left_join(rin, ck.ortho.cons)

# removed modified NTerm
rin <- rin %>% mutate(ortho_cons_ck = case_when( # REMOVE CCL5[5P7] MODIFIED NTERM
  file  == "5uiw" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5uiw" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ ortho_cons_ck,
  file  != "5uiw"  ~ ortho_cons_ck
))
rin <- rin %>% mutate(ortho_cons_ck = case_when( # REMOVE CX3CL1.35 MODIFIED NTERM
  file  == "5wb2" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5wb2" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ ortho_cons_ck,
  file  != "5wb2"  ~ ortho_cons_ck
))
rin <- rin %>% mutate(ortho_cons_ck = case_when( # REMOVE CCL5[6P4] MODIFIED NTERM
  file  == "7o7f" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "7o7f" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ ortho_cons_ck,
  file  != "7o7f"  ~ ortho_cons_ck
))

# combine protein-specific ortholog scores RECEPTOR
rin <- left_join(rin, ckr.ortho.cons)

# remove to clean
rm(ck.cons, ckr.cons, ck.ortho.cons, ckr.ortho.cons)


# (3) ADD CC/CXC PER PROTEIN PER POSITION SCORES -------------------------------
# add per protein, per protein cc/cxc classification scores CHEMOKINE
ck.cc.cxc.score <- read_csv("data/sequence/chemokine/processed/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  dplyr::select(protein, position, mean, sd) %>% unique()
ck.cc.cxc.score.virus <- read_csv("data/sequence/chemokine/processed/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  dplyr::select(protein, position, mean, sd) %>% unique() 
ck.cc.cxc.score.virus <- ck.cc.cxc.score.virus %>% filter(grepl("vmip2", ck.cc.cxc.score.virus$protein))
ck.cc.cxc.score.virus$protein <- c("vmipii")
ck.cc.cxc.score <- bind_rows(ck.cc.cxc.score, ck.cc.cxc.score.virus)
rm(ck.cc.cxc.score.virus)
colnames(ck.cc.cxc.score) <- c("ck", "source_gnccn", "cc_cxc_lr_score_ck", "cc_cxc_lr_score_sd_ck")
rin <- left_join(rin, ck.cc.cxc.score)
rm(ck.cc.cxc.score)

# add per protein, per protein cc/cxc classification scores RECEPTOR
ckr.cc.cxc.score <- read_csv("data/sequence/gpcr/processed/CKR_CLASSA_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  dplyr::select(protein, position, mean, sd) %>% unique()
ckr.cc.cxc.score.virus <- read_csv("data/sequence/gpcr/processed/CKR_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  dplyr::select(protein, position, mean, sd) %>% unique()
ckr.cc.cxc.score.virus <- ckr.cc.cxc.score.virus %>% filter(grepl("us28", ckr.cc.cxc.score.virus$protein))
ckr.cc.cxc.score.virus$protein <- c("us28")
ckr.cc.cxc.score <- bind_rows(ckr.cc.cxc.score, ckr.cc.cxc.score.virus)
rm(ckr.cc.cxc.score.virus)
colnames(ckr.cc.cxc.score) <- c("ckr", "target_gnccn", "cc_cxc_lr_score_ckr", "cc_cxc_lr_score_sd_ckr")
ckr.cc.cxc.score <- ckr.cc.cxc.score  %>% separate(target_gnccn, sep = "gn", into  = c("temp","target_gnccn"))
ckr.cc.cxc.score <- ckr.cc.cxc.score  %>% dplyr::select(-temp)
rin <- left_join(rin, ckr.cc.cxc.score)
rm(ckr.cc.cxc.score)

# removed modified NTerm
rin <- rin %>% mutate(cc_cxc_lr_score_ck = case_when( # REMOVE CCL5[5P7] MODIFIED NTERM
  file  == "5uiw" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5uiw" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cc_cxc_lr_score_ck,
  file  != "5uiw"  ~ cc_cxc_lr_score_ck
))
rin <- rin %>% mutate(cc_cxc_lr_score_ck = case_when( # REMOVE CX3CL1.35 MODIFIED NTERM
  file  == "5wb2" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5wb2" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cc_cxc_lr_score_ck,
  file  != "5wb2"  ~ cc_cxc_lr_score_ck
))
rin <- rin %>% mutate(cc_cxc_lr_score_ck = case_when( # REMOVE CCL5[6P4] MODIFIED NTERM
  file  == "7o7f" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "7o7f" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cc_cxc_lr_score_ck,
  file  != "7o7f"  ~ cc_cxc_lr_score_ck
))

rin <- rin %>% mutate(cc_cxc_lr_score_sd_ck = case_when( # REMOVE CCL5[5P7] MODIFIED NTERM
  file  == "5uiw" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5uiw" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cc_cxc_lr_score_sd_ck,
  file  != "5uiw"  ~ cc_cxc_lr_score_sd_ck
))
rin <- rin %>% mutate(cc_cxc_lr_score_sd_ck = case_when( # REMOVE CX3CL1.35 MODIFIED NTERM
  file  == "5wb2" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5wb2" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cc_cxc_lr_score_sd_ck,
  file  != "5wb2"  ~ cc_cxc_lr_score_sd_ck
))
rin <- rin %>% mutate(cc_cxc_lr_score_sd_ck = case_when( # REMOVE CCL5[6P4] MODIFIED NTERM
  file  == "7o7f" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "7o7f" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cc_cxc_lr_score_sd_ck,
  file  != "7o7f"  ~ cc_cxc_lr_score_sd_ck
))


# (4) ADD SNP AND CANCER -------------------------------------------------------
# import and add SNP/cancer chemokine
ck.snp.can <- read_csv("data/integrated/CK_CONS_CCCXC_SNP_CAN.csv") %>%
  dplyr::select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count) %>% unique()
colnames(ck.snp.can) <- c("ck", "source_gnccn", "dom1", "snp_count_ck", "snp_freq_count_ck", "cancer_mut_count_ck")
ck.snp.can$ck <- tolower(ck.snp.can$ck)
rin <- left_join(rin, ck.snp.can)
rm(ck.snp.can)

# import and add SNP/cancer receptor
ckr.snp.can <- read_csv("data/integrated/CKR_CONS_CCCXC_SNP_CAN.csv") %>%
  dplyr::select(protein, gn, dom, snp_count, snp_freq_count, cancer_mut_count) %>% unique()
colnames(ckr.snp.can) <- c("ckr", "target_gnccn", "dom2", "snp_count_ckr", "snp_freq_count_ckr", "cancer_mut_count_ckr")
ckr.snp.can$ckr <- tolower(ckr.snp.can$ckr)
ckr.snp.can <- ckr.snp.can  %>% separate(target_gnccn, sep = "gn", into  = c("temp","target_gnccn"))
ckr.snp.can <- ckr.snp.can  %>% dplyr::select(-temp)
rin <- left_join(rin, ckr.snp.can)
rm(ckr.snp.can)

# removed modified Nterm
rin <- rin %>% mutate(snp_count_ck = case_when( # REMOVE CCL5[5P7] MODIFIED NTERM
  file  == "5uiw" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5uiw" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ snp_count_ck,
  file  != "5uiw"  ~ snp_count_ck
))
rin <- rin %>% mutate(snp_count_ck = case_when( # REMOVE CX3CL1.35 MODIFIED NTERM
  file  == "5wb2" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5wb2" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ snp_count_ck,
  file  != "5wb2"  ~ snp_count_ck
))
rin <- rin %>% mutate(snp_count_ck = case_when( # REMOVE CCL5[6P4] MODIFIED NTERM
  file  == "7o7f" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "7o7f" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ snp_count_ck,
  file  != "7o7f"  ~ snp_count_ck
))
rin <- rin %>% mutate(snp_freq_count_ck = case_when( # REMOVE CCL5[5P7] MODIFIED NTERM
  file  == "5uiw" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5uiw" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ snp_freq_count_ck,
  file  != "5uiw"  ~ snp_freq_count_ck
))
rin <- rin %>% mutate(snp_freq_count_ck = case_when( # REMOVE CX3CL1.35 MODIFIED NTERM
  file  == "5wb2" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5wb2" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ snp_freq_count_ck,
  file  != "5wb2"  ~ snp_freq_count_ck
))
rin <- rin %>% mutate(snp_freq_count_ck = case_when( # REMOVE CCL5[6P4] MODIFIED NTERM
  file  == "7o7f" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "7o7f" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ snp_freq_count_ck,
  file  != "7o7f"  ~ snp_freq_count_ck
))
rin <- rin %>% mutate(cancer_mut_count_ck = case_when( # REMOVE CCL5[5P7] MODIFIED NTERM
  file  == "5uiw" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5uiw" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cancer_mut_count_ck,
  file  != "5uiw"  ~ cancer_mut_count_ck
))
rin <- rin %>% mutate(cancer_mut_count_ck = case_when( # REMOVE CX3CL1.35 MODIFIED NTERM
  file  == "5wb2" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "5wb2" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cancer_mut_count_ck,
  file  != "5wb2"  ~ cancer_mut_count_ck
))
rin <- rin %>% mutate(cancer_mut_count_ck = case_when( # REMOVE CCL5[6P4] MODIFIED NTERM
  file  == "7o7f" & grepl("NTc.Cm", rin$source_gnccn) ~  NA_real_,
  file  == "7o7f" & (!grepl("NTc.Cm", rin$source_gnccn)) ~ cancer_mut_count_ck,
  file  != "7o7f"  ~ cancer_mut_count_ck
))

# add MARKING columns
rin$STRUCTURE <- c("STRUCTURE")
rin$CONSERVATION <- c("CONSERVATION")
rin$CLASSIFICATION <- c("CLASSIFICATION")
rin$POLYMORPHISM <- c("POLYMORPHISM")


# reorder to pair CK and CKR  scores
rin <- rin %>%
  dplyr::select(file, ck, ckr, source_gnccn, target_gnccn,  dom1,  dom2, 
         STRUCTURE,
         no_pdb,
         no_contacts_file_ck,  no_contacts_file_ckr, no_contacts_all_ck, no_contacts_all_ckr,
         CONSERVATION,
         all_para_ck, all_para_ckr, all_non_ackr_para_ckr,
         all_classa_ckr,
         all_cc_cxc_para_ck, all_cc_cxc_para_ckr, 
         cc_para_ck, cc_para_ckr,
         cxc_para_ck, cxc_para_ckr,
         ack_para_ckr,
         ortho_cons_ck, ortho_cons_ckr,
         CLASSIFICATION,
         cc_cxc_lr_ck, cc_cxc_lr_ckr,
         cc_cxc_lr_sd_ck, cc_cxc_lr_sd_ckr,
         cc_cxc_lr_score_ck, cc_cxc_lr_score_ckr,
         cc_cxc_lr_score_sd_ck, cc_cxc_lr_score_sd_ckr,
         POLYMORPHISM,
         snp_count_ck, snp_count_ckr,
         snp_freq_count_ck, snp_freq_count_ckr,
         cancer_mut_count_ck, cancer_mut_count_ckr)

# write
# write_csv(rin, "data/integrated/RIN_CONS_CLASS.csv") # LAST WRITTEN 20231226

