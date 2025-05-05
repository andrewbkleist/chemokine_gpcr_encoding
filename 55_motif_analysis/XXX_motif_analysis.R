# Name:     02_motif_analysis.R
# Updated:  20210510
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)

##### XXX ######################################################################

  ck <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")
  ckr.nt <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")
  ckr.ecl2 <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_ECL2.csv")

  
  
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") 
  rin <- rin %>% filter(file == "ngo") 
  satmut <- read_csv("09_sat_mut_data/output/cxcr4_clean_gpcrdb_means_log.csv") %>% 
    select(gn, sele, sub_mean) %>%
    filter(sele == "cxcl12") %>%
    group_by(gn) %>%
    summarise(sub_mean2 = mean(sub_mean)) %>%
    ungroup()
  colnames(satmut)[1] <- c("target_gnccn")
  rin <- left_join(rin, satmut)
  
  
  cxcl12 <- rin %>%
    filter(dom1 == "NTc" | dom2 == "NTr") %>%
    filter(file == "ngo") %>%
    ggplot(aes(all_non_ackr_para_ckr, all_para_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
    # scale_size(range = c(2, 8)) +
    xlim(0,1) +
    ylim(0,1) +
    theme_minimal()
  
  
  rin <- rin %>%
    filter(dom1 == "NTc" | dom2 == "NTr") %>%
    filter(file == "6wwz")
  
    ggplot(aes(all_non_ackr_para_ckr, all_para_ck))  +
    geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
    # scale_size(range = c(2, 8)) +
    xlim(0,1) +
    ylim(0,1) +
    theme_minimal()
  
  
  
#### XXX ######################################################################
  
  # correlation with promiscuity...
  ck <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")
  ck <- ck %>% filter(protein %in%  c("ccl27", "ccl28", "ccl2", "ccl5", "ccl7", "ccl8", "ccl11", "ccl13", "ccl14", "ccl15"))
  ckr.nt <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_NTERM.csv")
  ckr.ecl2 <- read_csv("08_ckr_motif/output/CKR_MOTIF_FREQUENCY_ECL2.csv")
  
  # read how promiscuous
  prom <- read_csv("06_network/output/chemokine_gpcr_interactome_validated.csv")
  ck.prom <- prom %>% count(ck)  
  ckr.prom <- prom %>% count(ckr)
  colnames(ck.prom)[1] <- c("protein")
  colnames(ckr.prom)[1] <- c("protein")
  
  ck.prom$protein <- tolower(ck.prom$protein)
  ckr.prom$protein <- tolower(ckr.prom$protein)
  
  # map
  ck <- left_join(ck, ck.prom)
  ckr.nt <- left_join(ckr.nt, ckr.prom)
  ckr.ecl2 <- left_join(ckr.ecl2, ckr.prom)
  
  
  ckr.ecl2 %>%
    filter(mer == "mer2") %>%
    ggplot(aes(as.character(n), pct_super)) +
    geom_boxplot() +
    theme_minimal()
  
  
  
  