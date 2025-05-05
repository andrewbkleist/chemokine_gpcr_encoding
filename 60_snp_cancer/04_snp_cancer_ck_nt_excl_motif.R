# Name:     01_map_motif_to_seq.R
# Updated:  20210703
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

# NOTE 1 (20210703): ASSOCIATING MOTIFS WITH SPECIFIC POSITIONS
# This script will associate motifs with PER POSITION cancer and SNP counts. 
# This is complicated by the fact that motifs comprise multiple residues whereas
# SNPs and cancer are associated with individual residue positions. The question
# is the following: how can you associate SNPs and mutations (individual residues)
# strings of multiple residues? Two ways are the following (1) keep the dataframe
# with respect to individual residues, and then label the residue position
# by whether it belongs to a motif or not, or (2) convert canncer/SNP to new
# dataframe where motif is the first column and canncer/SNP values associated with
# the residues of that motif are listed (e.g. avergae no. SNPs for all residues
# compirising that motif). Approach (1) is more straightforward but must address
# how to define residues that are part of both "motifs" (>50% cons.) and "non-motifs"
# (<50% cons).  To address, will distinguish "exclusive motif" residues,
# which are residues that ONLY belong to motifs. Did at 3-mer level.


##### 0: IDENTIFY EXCLUSIVE NTERM ##############################################

  # (0.1) ASSOCIATE MOTIFS WITH ALN POSITION NOMENCLATURE ----------------------
  #   (see /6_snp_cancer/3_get_motifs.R for original script)

  # import MOTIFS with mapped positional information
  lookup.motif <- read_csv("07_ck_motif/data/processed/20210708_motif_conversion_3mer_ck_nterm.csv")
  
  # import MOTIF CONSERVATION
  motif.cons <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv") %>% 
    #filter(pct_ortho >= 0.5) %>% 
    filter(mer == "mer3")
  motif.cons <- motif.cons %>% dplyr::select(motif, protein, pct_ortho)
  motif.cons$protein <- toupper(motif.cons$protein)
  
  # JOIN ...
  lookup.motif <- left_join(lookup.motif, motif.cons)
  lookup.motif <- lookup.motif %>% filter(!is.na(motif))
  lookup.motif <- lookup.motif %>% mutate(tier = case_when(
    pct_ortho >= 0.5 ~ "motif",
    pct_ortho < 0.5 ~ "fragment"
  ))
  rm(motif.cons)
  
  # (0.2) ISOLATE EXCLUSIVE MOTIF POSITIONS ------------------------------------
  lookup.motif <- lookup.motif %>%
    dplyr::select(-motif, -pct_ortho) %>% unique()
  lookup.motif <- lookup.motif %>% mutate(value = case_when(
    tier ==  "motif" ~ 1,
    tier ==  "fragment" ~ 1
  ))
  lookup.motif <- lookup.motif %>% filter(!is.na(tier)) # removes XCL1, XCL2, CCL4L1
  lookup.motif <- lookup.motif %>% spread(tier, value)
  lookup.motif[is.na(lookup.motif)] <- 0
  
  lookup.motif <- lookup.motif %>% mutate(frag_or_no = case_when(
    fragment == 0 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 0 ~ "fragment",
    fragment == 1 & motif == 1 ~ "fragment",
    fragment == 0 & motif == 1 ~ "motif"
  ))
  lookup.motif <- lookup.motif %>% dplyr::select(-motif, -fragment)
  
  write_csv(lookup.motif, "60_snp_cancer/data/processed/ck_excl_motif_aln_pos.csv") 
  
  
##### 1: CHEMOKINE NTERM SNP ######################################################
  
  # (1.1) IMPORT, COMBINE, CHI SQUARE ------------------------------------------
  # import lookup 
  lookup.motif <- read_csv("60_snp_cancer/data/processed/ck_excl_motif_aln_pos.csv")
  
  # import, reformat
  data <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")
  
  # filter extreme N-term, batch approach (ie remove N-termini >10 residues)
  # & remove C-terminus
  data <- data %>% filter(!(ccn %in% c(
    "NTc.Cm11","NTc.Cm12","NTc.Cm13","NTc.Cm14","NTc.Cm15","NTc.Cm16","NTc.Cm17","NTc.Cm18","NTc.Cm19","NTc.Cm20","NTc.Cm21","NTc.Cm22","NTc.Cm23","NTc.Cm24","NTc.Cm25","NTc.Cm26","NTc.Cm27","NTc.Cm28","NTc.Cm29","NTc.Cm30","NTc.Cm31","NTc.Cm32","NTc.Cm33","NTc.Cm34","NTc.Cm35","NTc.Cm36","NTc.Cm37","NTc.Cm38","NTc.Cm39","NTc.Cm40","NTc.Cm41","NTc.Cm42","NTc.Cm43","NTc.Cm44","NTc.Cm45","NTc.Cm46","NTc.Cm47","NTc.Cm48","NTc.Cm49","NTc.Cm50","NTc.Cm51","NTc.Cm52","NTc.Cm53","NTc.Cm54","NTc.Cm55","NTc.Cm56","NTc.Cm57","NTc.Cm58","NTc.Cm59","NTc.Cm60","NTc.Cm61","NTc.Cm62","NTc.Cm63","NTc.Cm64","NTc.Cm65","NTc.Cm66","NTc.Cm67","NTc.Cm68","NTc.Cm69","NTc.Cm70","NTc.Cm71","NTc.Cm72","NTc.Cm73","NTc.Cm74","NTc.Cm75","NTc.Cm76","NTc.Cm77","NTc.Cm78","NTc.Cm79","NTc.Cm80","NTc.Cm81","NTc.Cm82","NTc.Cm83","NTc.Cm84","NTc.Cm85","NTc.Cm86","NTc.Cm87","NTc.Cm88","NTc.Cm89","NTc.Cm90","NTc.Cm91","NTc.Cm92","NTc.Cm93","NTc.Cm94","NTc.Cm95"
  ))) %>% filter(dom != "CT")
  
  
  # select relevant cols
  data <- data %>% select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count)

  # filter high frequency alleles (30,000, ~10% frequency in population)
  data <- data %>% filter(snp_freq_count < 30000)
  
  # map "motif_or_no "exclusive motif" positions to new df
  data <- left_join(data, lookup.motif)
  
  # get numbers of snps in exlcusive motif positions versus "fragment" positions
  nt <- data %>% 
    filter(!is.na(frag_or_no)) %>%
    group_by(frag_or_no) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(snp_freq_count)) %>% ungroup()
  chisq.test(c(3633, 13744), p=c(164, 206)/(164 + 206))
  
  
  # (1.2) GRAPH RESULTS --------------------------------------------------------

  # define function to graph expected and observed chi square values for both groups
  ChiToGraph <- function(DF, TIER, NON){
    data <- DF %>% gather(key, value, 2:3)
    data <- data %>% dplyr::group_by(key) %>% dplyr::mutate(total = sum(value)) %>% dplyr::ungroup()
    data <- data %>% mutate(fraction = value / total)
    data <- data.frame("EO" = c("O", "E", "O", "E"), 
                       "tier" = c(TIER, TIER, NON, NON),
                       "value" = c(data$value[4], data$fraction[2]*data$total[3], data$value[3], data$fraction[1]*data$total[3]) )
    return(data)
    rm(data)
  }
  
  # PLOT - EXPECTED VS. OBSERVED
  test <- ChiToGraph(nt, "motif", "fragment")
  
  order = c("motif", "fragment")
  order2 = c("E", "O")
  test$tier <- factor(test$tier, levels = order)
  test$EO <- factor(test$EO, levels = order2)
  test %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # Hypothesis: chemokines have more SNPs than expected in N-terminus
  # Results: True! When you remove C-terminal region
  
  # (best would be to then show unique activity for Nterm SNP/cancer)
  # coulds show variation in Nterm (ideally motif, ie conserved) leads to
  # functional innovation in human DISEASE
  
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    ylim(-5, 0.6)
    
  
  
##### 2: CHEMOKINE NTERM CANCER ################################################
  
  # (2.1) IMPORT, COMBINE, CHI SQUARE ------------------------------------------
  # import lookup 
  lookup.motif <- read_csv("60_snp_cancer/data/processed/ck_excl_motif_aln_pos.csv")
  # import, reformat
  data <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")
  
  # filter extreme N-term, batch approach (ie remove N-termini >10 residues)
  # & remove C-terminus
  data <- data %>% filter(!(ccn %in% c(
    "NTc.Cm11","NTc.Cm12","NTc.Cm13","NTc.Cm14","NTc.Cm15","NTc.Cm16","NTc.Cm17","NTc.Cm18","NTc.Cm19","NTc.Cm20","NTc.Cm21","NTc.Cm22","NTc.Cm23","NTc.Cm24","NTc.Cm25","NTc.Cm26","NTc.Cm27","NTc.Cm28","NTc.Cm29","NTc.Cm30","NTc.Cm31","NTc.Cm32","NTc.Cm33","NTc.Cm34","NTc.Cm35","NTc.Cm36","NTc.Cm37","NTc.Cm38","NTc.Cm39","NTc.Cm40","NTc.Cm41","NTc.Cm42","NTc.Cm43","NTc.Cm44","NTc.Cm45","NTc.Cm46","NTc.Cm47","NTc.Cm48","NTc.Cm49","NTc.Cm50","NTc.Cm51","NTc.Cm52","NTc.Cm53","NTc.Cm54","NTc.Cm55","NTc.Cm56","NTc.Cm57","NTc.Cm58","NTc.Cm59","NTc.Cm60","NTc.Cm61","NTc.Cm62","NTc.Cm63","NTc.Cm64","NTc.Cm65","NTc.Cm66","NTc.Cm67","NTc.Cm68","NTc.Cm69","NTc.Cm70","NTc.Cm71","NTc.Cm72","NTc.Cm73","NTc.Cm74","NTc.Cm75","NTc.Cm76","NTc.Cm77","NTc.Cm78","NTc.Cm79","NTc.Cm80","NTc.Cm81","NTc.Cm82","NTc.Cm83","NTc.Cm84","NTc.Cm85","NTc.Cm86","NTc.Cm87","NTc.Cm88","NTc.Cm89","NTc.Cm90","NTc.Cm91","NTc.Cm92","NTc.Cm93","NTc.Cm94","NTc.Cm95"
  ))) %>% filter(dom != "CT")
  
  # select relevant cols
  data <- data %>% select(protein, ccn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # filter high frequency alleles (30,000, ~10% frequency in population)
  # data <- data %>% filter(snp_freq_count < 30000)
  
  # map "motif_or_no "exclusive motif" positions to new df
  data <- left_join(data, lookup.motif)
  
  # get numbers of snps in exlcusive motif positions versus "fragment" positions
  nt <- data %>% 
    filter(!is.na(frag_or_no)) %>%
    group_by(frag_or_no) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(cancer_mut_count)) %>% ungroup()
  chisq.test(c(24, 54), p=c(188, 318)/(188 + 318))
  
  
  # (2.2) GRAPH RESULTS --------------------------------------------------------
  
  # define function to graph expected and observed chi square values for both groups
  ChiToGraph <- function(DF, TIER, NON){
    data <- DF %>% gather(key, value, 2:3)
    data <- data %>% dplyr::group_by(key) %>% dplyr::mutate(total = sum(value)) %>% dplyr::ungroup()
    data <- data %>% mutate(fraction = value / total)
    data <- data.frame("EO" = c("O", "E", "O", "E"), 
                       "tier" = c(TIER, TIER, NON, NON),
                       "value" = c(data$value[4], data$fraction[2]*data$total[3], data$value[3], data$fraction[1]*data$total[3]) )
    return(data)
    rm(data)
  }
  
  # PLOT - EXPECTED VS. OBSERVED
  test <- ChiToGraph(nt, "motif", "fragment")
  
  order = c("motif", "fragment")
  order2 = c("E", "O")
  test$tier <- factor(test$tier, levels = order)
  test$EO <- factor(test$EO, levels = order2)
  test %>%
    ggplot(aes(tier, value, group=EO, fill = tier)) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_manual(values=c("grey30", "grey70")) +
    theme_minimal()
  
  # Hypothesis: chemokines have more SNPs than expected in N-terminus
  # Results: True! When you remove C-terminal region
  
  # (best would be to then show unique activity for Nterm SNP/cancer)
  # coulds show variation in Nterm (ideally motif, ie conserved) leads to
  # functional innovation in human DISEASE
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() 
  