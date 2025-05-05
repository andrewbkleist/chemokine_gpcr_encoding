# Name:     02_integrate_rins.R
# Updated:  20220929
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE 1 (20220930): Choosing complexes for consensus contact prep
#
# Summary of complexes (grouped by type/receptor):
# (1) CCL5 : CCR5 & VARIANTS (SEMI-DEGENERATE)
# 5uiw: ccl5-5p7 / ccr5
# 7o7f: ccl5-6p4 / ccr5 **
# 7f1r: ccl5 / ccr5 (incomplete) **
# zheng: ccl5 / ccr5 (model)
# 
# (2) OTHER CCL* : CCR* COMPLEXES
# 7vl9: ccl15 / ccr1
# 7xa3: ccl2 / ccr2
# 7f1t: ccl3 / ccr5 
# 6wwz: ccl20 / ccr6
#
# (3) CXCL* : CXCR* COMPLEXES
# 6lfo: cxcl8 / cxcr2
# ngo: cxcl12 / cxcr4  (model)
#
# (4) NON-CANNONICAL HUMAN COMPLEXES
# 7sk3: cxcl12 / ackr3
# 7xbx: cx3cl1 / cx3cr1
#
# (5) VIRAL CK:CKR COMPONENTS
# 4rws: vmipii / cxcr4
# 4xt1: cx3cl1 / us28
# 5wb2: cx3cl1.35 / us28 **
#
# (6) NON- CK/CKR COMPONENTS
# 6meo: gp120 / ccr5 (HIV complex)
# ...could also add receptor-drug complexes
#
# Total: 15 complexes (not including 6meo)
# Non-degenerate: 11 complexes
# Non-degenerate human: 9 complexes
#
# 
# OPTIONS:
# (1) Include ALL complexes
# PROS: Semi-objective since its not curated; CONS: degenerate complexes
# includiing 4x CCL5-CCR5 and 2x CX3CL1:US28. OTHER: Could do second version in 
# supp in which degenerates are removed to evaluate that "consensus" contacts are
# not artificially enriched because of degeneracy. Would still remove HIV-CCR5.
# Would also not want to make explicit threshold for what is considered
# consensus versus not, could evaluate each abundant position independently
# (eg CCL5:CCR5 native complex is incomplete). Another "pro" is that the 
# "degenerate" complexes are not truly degenerate because of distinct N-termini
# and thus sample unique chemokine-GPCR contacts "allowed" by the architecture.
#
# (2) Include NON-DEGENERATE complexes
# PROS: Only represents each complex 1x (eg choose 1x CCL5-CCR5 and 1x CX3CL1-
# US28). Top contacts would be more likely representative of  family. CONS: 
# Would have to choose among 4x CCL5:CCR5 complexes, all of which are individually
# non-ideal (and actually distinct in N-terminal regions, see above).
#
# (3) MULTIPLE overlapping sets
# Variation of (1) and (2) in which you include non-degenerate complexes and
# swap in different sets of the redundant complexes and create some sort of
# summary statistic that is representative of the various sets (eg mean no.
# contacts per position).
#
# Below code persues option (1) above

    
##### 1: COMBINE PAIRWISE CONTACTS AND CONSERVATION - 0.5 ANGSTROM (DEFAULT) ###
  
  # (2.1) Contact analysis  ----------------------------------------------------
    # import rin
    rin <- read_csv("01_structure_contacts/output/RIN_residue.csv") %>%
      filter(class == "full") %>% filter(Chain1 != Chain2) %>%
      select(file, source_gnccn, target_gnccn, dom1, dom2)
    
    # add representation of  each contact among *all* complexes
    n.pdb <- rin %>% select(source_gnccn, target_gnccn) %>% count(source_gnccn, target_gnccn)
    colnames(n.pdb)[3] <- c("no_pdb")
    rin <- left_join(rin, n.pdb)
    rm(n.pdb)
    
    # add number of unique residues contacted by each chemokine/GPCR residue WITHIN COMPLEX
    n.resi.ck <- rin %>% select(file, source_gnccn) %>% count(file, source_gnccn)
    colnames(n.resi.ck)[3] <- c("no_contacts_file_ck")
    n.resi.ckr <- rin %>% select(file, target_gnccn) %>% count(file, target_gnccn)
    colnames(n.resi.ckr)[3] <- c("no_contacts_file_ckr")
    rin <- left_join(rin, n.resi.ck)
    rin <- left_join(rin, n.resi.ckr)
    
    rm(n.resi.ck, n.resi.ckr)
    
    # add number of unique residues contacted by each chemokine/GPCR residue ACROSS COMPLEXES
    n.resi.ck <- rin %>% 
      filter(file != "6meo") %>%
      # remove gp120-CCR5 complex so as not to skew contacts
      select(source_gnccn, target_gnccn) %>% unique() %>% 
      select(source_gnccn) %>% count(source_gnccn)
    colnames(n.resi.ck)[2] <- c("no_contacts_all_ck")
    
    n.resi.ckr <- rin %>% 
      filter(file != "6meo") %>%
      # remove gp120-CCR5 complex so as not to skew contacts
      select(source_gnccn, target_gnccn) %>% unique() %>% 
      select(target_gnccn) %>% count(target_gnccn) 
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
      file == "7xa3" ~ "ccl2"
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
      file == "7xa3" ~ "ccr2"
    ))
    rin <- rin %>% select(file, ck, ckr, source_gnccn, target_gnccn, dom1, dom2, 
                  no_pdb, no_contacts_file_ck, no_contacts_file_ckr,
                  no_contacts_all_ck, no_contacts_all_ckr)
    
    rm(n.resi.ck, n.resi.ckr)
  
    
  # (2.2) ADD CONSERVATION -----------------------------------------------------
    
    # import conservation chemokine
    ck.cons  <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")
    colnames(ck.cons) <- c("ck", "source_gnccn","dom1", "all_para_ck", "all_cc_cxc_para_ck",
                           "cc_para_ck", "cxc_para_ck", "ortho_cons_ck",  "ngaps_ck",  "pctgaps_ck",
                           "cc_cxc_lr_ck", "cc_cxc_lr_sd_ck", "cc_cxc_lr_score_ck", "cc_cxc_lr_score_sd_ck",
                           "snp_count_ck", "snp_freq_count_ck", "cancer_mut_count_ck")
    ck.ortho.cons <- ck.cons %>% select(ck, source_gnccn, ortho_cons_ck) %>% 
      separate(ck, sep = "_", into = "ck")
    ck.cons <- ck.cons %>% 
      select("source_gnccn","dom1", "all_para_ck",  "all_cc_cxc_para_ck", 
             "cc_para_ck", "cxc_para_ck", "cc_cxc_lr_ck", "cc_cxc_lr_sd_ck") %>%
      unique()
    
    # import conservation receptor
    ckr.cons  <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv")  
    colnames(ckr.cons) <- c("ckr", "target_gnccn","dom2", "all_para_ckr", "all_non_ackr_para_ckr", "all_cc_cxc_para_ckr", "all_classa_ckr",  
                           "cc_para_ckr", "cxc_para_ckr", "ack_para_ckr", "ortho_cons_ckr",  "ngaps_ckr",  "pctgaps_ckr",
                           "cc_cxc_lr_ckr", "cc_cxc_lr_sd_ckr", "cc_cxc_lr_score_ckr", 
                           "cc_cxc_lr_score_sd_ckr",
                           "snp_count_ckr", "snp_freq_count_ckr", 
                           "cancer_mut_count_ckr")
    ckr.cons <- ckr.cons  %>% separate(target_gnccn, sep = "gn", into  = c("temp","target_gnccn"))
    ckr.cons <- ckr.cons  %>% select(-temp)
    
    ckr.ortho.cons <- ckr.cons %>% select(ckr, target_gnccn, ortho_cons_ckr) %>% 
      separate(ckr, sep = "_", into = "ckr")
    ckr.cons <- ckr.cons %>% 
      select("target_gnccn","dom2", "all_para_ckr", "all_non_ackr_para_ckr", "all_cc_cxc_para_ckr","all_classa_ckr",  
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
    
    
  # (2.3) ADD CC/CXC PER PROTEIN PER POSITION SCORES ---------------------------
    
    # add per protein, per protein cc/cxc classification scores CHEMOKINE
    ck.cc.cxc.score <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      select(protein, position, mean, sd) %>% unique()
    ck.cc.cxc.score.virus <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      select(protein, position, mean, sd) %>% unique() 
    ck.cc.cxc.score.virus <- ck.cc.cxc.score.virus %>% filter(grepl("vmip2", ck.cc.cxc.score.virus$protein))
    ck.cc.cxc.score.virus$protein <- c("vmipii")
    ck.cc.cxc.score <- bind_rows(ck.cc.cxc.score, ck.cc.cxc.score.virus)
    rm(ck.cc.cxc.score.virus)
    colnames(ck.cc.cxc.score) <- c("ck", "source_gnccn", "cc_cxc_lr_score_ck", "cc_cxc_lr_score_sd_ck")
    rin <- left_join(rin, ck.cc.cxc.score)
    rm(ck.cc.cxc.score)
    
    # add per protein, per protein cc/cxc classification scores RECEPTOR
    ckr.cc.cxc.score <- read_csv("03_ckr_seq/output/CKR_CLASSA_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      select(protein, position, mean, sd) %>% unique()
    ckr.cc.cxc.score.virus <- read_csv("03_ckr_seq/output/CKR_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      select(protein, position, mean, sd) %>% unique()
    ckr.cc.cxc.score.virus <- ckr.cc.cxc.score.virus %>% filter(grepl("us28", ckr.cc.cxc.score.virus$protein))
    ckr.cc.cxc.score.virus$protein <- c("us28")
    ckr.cc.cxc.score <- bind_rows(ckr.cc.cxc.score, ckr.cc.cxc.score.virus)
    rm(ckr.cc.cxc.score.virus)
    colnames(ckr.cc.cxc.score) <- c("ckr", "target_gnccn", "cc_cxc_lr_score_ckr", "cc_cxc_lr_score_sd_ckr")
    ckr.cc.cxc.score <- ckr.cc.cxc.score  %>% separate(target_gnccn, sep = "gn", into  = c("temp","target_gnccn"))
    ckr.cc.cxc.score <- ckr.cc.cxc.score  %>% select(-temp)
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
    
    
  # (2.4) ADD SNP AND CANCER ---------------------------------------------------
    # import and add SNP/cancer chemokine
    ck.snp.can <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv") %>%
      select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count) %>% unique()
    colnames(ck.snp.can) <- c("ck", "source_gnccn", "dom1", "snp_count_ck", "snp_freq_count_ck", "cancer_mut_count_ck")
    ck.snp.can$ck <- tolower(ck.snp.can$ck)
    rin <- left_join(rin, ck.snp.can)
    rm(ck.snp.can)
    
    # import and add SNP/cancer receptor
    ckr.snp.can <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv") %>%
      select(protein, gn, dom, snp_count, snp_freq_count, cancer_mut_count) %>% unique()
    colnames(ckr.snp.can) <- c("ckr", "target_gnccn", "dom2", "snp_count_ckr", "snp_freq_count_ckr", "cancer_mut_count_ckr")
    ckr.snp.can$ckr <- tolower(ckr.snp.can$ckr)
    ckr.snp.can <- ckr.snp.can  %>% separate(target_gnccn, sep = "gn", into  = c("temp","target_gnccn"))
    ckr.snp.can <- ckr.snp.can  %>% select(-temp)
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
      select(file, ck, ckr, source_gnccn, target_gnccn,  dom1,  dom2, 
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
  # write_csv(rin, "05_integrate/output/RIN_CONS_CLASS.csv") # LAST WRITTEN 20221004
  