# Name:     05_snp_cancer_ckr_ecl2_excl_motif.R
# Updated:  20210712
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### 0: IDENTIFY EXCLUSIVE ECL2 ###############################################
  
  # (0.1) ASSOCIATE MOTIFS WITH ALN POSITION NOMENCLATURE ----------------------
  #   (see /6_snp_cancer/3_get_motifs.R for original script)
  
  # import MOTIFS with mapped positional information
  lookup.motif <- read_csv("08_ckr_motif/data/processed/20210710_motif_conversion_3mer_ckr_ecl2.csv")
  
  # import MOTIF CONSERVATION
  motif.cons <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_ECL2.csv") %>% 
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
  
  write_csv(lookup.motif, "60_snp_cancer/data/processed/ckr_excl_motif_ecl2_aln_pos.csv") 


##### 1: RECEPTOR ECL2 SNP ######################################################
  
  # (1.1) IMPORT, COMBINE, CHI SQUARE ------------------------------------------
  # import lookup 
  lookup.motif <- read_csv("60_snp_cancer/data/processed/ckr_excl_motif_ecl2_aln_pos.csv")
  
  # import, reformat
  data <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv")
  
  # select relevant cols
  data <- data %>% select(protein, gn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # filter high frequency alleles (30,000, ~10% frequency in population)
  data <- data %>% filter(snp_freq_count < 30000)
  
  # map "motif_or_no "exclusive motif" positions to new df
  data <- left_join(data, lookup.motif)
  
  # get numbers of snps in exlcusive motif positions versus "fragment" positions
  nt <- data %>% 
    filter(!is.na(frag_or_no)) %>%
    group_by(frag_or_no) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(snp_freq_count)) %>% ungroup()
  chisq.test(c(367, 678), p=c(101, 237)/(101 + 237))
  
  
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
  
  # PLOT - EXPECTED VS. OBSERVED LOG RATIO PLOT
  test <- test %>% pivot_wider(names_from = EO, values_from = value)
  test <- test %>% mutate(oe_ratio = log2(O/E))
  test %>%
    ggplot(aes(tier, oe_ratio)) +
    geom_bar(stat = "identity") +
    theme_minimal() 


##### 2: RECEPTOR ECL2 CANCER ##################################################
  
  # (2.1) IMPORT, COMBINE, CHI SQUARE ------------------------------------------
  # import lookup 
  lookup.motif <- read_csv("60_snp_cancer/data/processed/ckr_excl_motif_ecl2_aln_pos.csv")
  
  # import, reformat
  data <- read_csv("05_integrate/output/CKR_CONS_CCCXC_SNP_CAN.csv")
  
  # select relevant cols
  data <- data %>% select(protein, gn, dom, snp_count, snp_freq_count, cancer_mut_count)
  
  # filter high frequency alleles (30,000, ~10% frequency in population)
  # data <- data %>% filter(snp_freq_count < 30000)
  
  # map "motif_or_no "exclusive motif" positions to new df
  data <- left_join(data, lookup.motif)
  
  # get numbers of snps in exlcusive motif positions versus "fragment" positions
  nt <- data %>% 
    filter(!is.na(frag_or_no)) %>%
    group_by(frag_or_no) %>%
    dplyr::summarise(n=dplyr::n(), count = sum(cancer_mut_count)) %>% ungroup()
  chisq.test(c(25, 62), p=c(101, 237)/(101 + 237))
  
  
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

  