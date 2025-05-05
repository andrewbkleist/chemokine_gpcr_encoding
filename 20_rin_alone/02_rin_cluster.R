# Name:     02_rin_cluster.R
# Updated:  20221007
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)


################################################################################

  # (1.1) COUNT CONTACTS  ------------------------------------------------------
  # import contacts, select interface only, select Xray only
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # (1.2) MAP TO MATRIX, HCLUST, DENDROGRAM --------------------------------------------
  rin <- rin %>% select(file, source_gnccn, target_gnccn) 
  
  rin <- rin %>% unite(rin, c(source_gnccn, target_gnccn), sep = "_")
  rin$n <- c(1)

  # spread
  rin <- rin %>% pivot_wider(names_from = rin, values_from = n)
  
  # gather, then replace NULL with "0"
  rin <- rin %>% pivot_longer(cols = c("NTc.Cm10_2x60":"b1b2.16_5x33"), names_to = "rin", values_to =  "n")
  rin$n <- as.character((rin$n))
  rin <- rin %>% mutate(n = case_when(
    n == "NULL" ~ "0",
    n == "1" ~ n
  ))
  rin$n <- as.numeric((rin$n))
  
  # select complexes (removing 6meo and 7f1r)
  complexes = c("5uiw", "4rws", "4xt1", "5wb2", "6lfo", "6wwz",
                "ngo", "zheng", "7xbx", "7vl9", "7sk3", "7o7f",
                "7f1t", "7xa3")
  
  # complexes = c("5uiw", "7o7f", "zheng", "7vl9", "7xa3", "7f1t",
  #               "6wwz", "6lfo", "ngo")
  
  rin <- rin %>% filter(file %in% complexes)
  rm(complexes)
  
  # re-spread
  rin <- rin %>% pivot_wider(names_from = rin, values_from = n)
  
  # spread again
  rownames(rin) <- rin$file
  rin <- as.matrix(rin)
  rin <- rin[,2:ncol(rin)]
  
  # similarity
  rin.dist <- dist(rin, method = "euclidean")
  hc <- hclust(rin.dist, method = "ward.D2")
  plot(hc,  hang = -1, cex = 0.6)
  
    # SUMMARY
    # CLUSTER 1 (CCR5, CCR1)
    # 7vl9: ccl15:ccr1; 7o7f: ccl5-6p4:ccr5; 7f1t: ccl3:ccr5; 5uiw: ccl5-5p7:ccr5
    # zheng: ccl5 / ccr5 (model)
    
    # CLUSTER 2
    # ngo: cxcl12:cxcr4  (model)
    
    # CLUSTER 3 (CCR2, CCR6)
    # 7xa3: ccl2:ccr2; 4rws: vmipii:cxcr4; 6wwz: ccl20:ccr6
    
    # CLUSTER 4
    # 7sk3: cxcl12:ackr3
    
    # CLUSTER 5 (US28, CX3CR1, CXCR2)
    # 4xt1: cx3cl1:us28; 7xbx: cx3cl1:cx3cr1; 5wb2: cx3cl1.35:us28; 6lfo: cxcl8:cxcr2
    
    # Generalized clustering of contacts made by ligand-receptor complexes 
    # within networks, consistent with residue-residue contacts driving
    # coupling networks; result functions similarly to conservation diversity despite
    # structural similarity results, provides rationale for investigating structure
    # to help explain selectivity. Also some surprises where diverse N-terminal
    # contacts seem to "rewire" coupling networks away from those of related
    # chemokines/receptors (e.g. CCL5-6P4 is clustered with CCL15-CCR1 complex)
  

  