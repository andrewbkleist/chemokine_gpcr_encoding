source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

temp <- as.data.frame(NULL)
seq <- c("ccr1","ccr2","ccr3","ccr4","ccr5","ccr6","ccr7","ccr8","ccr9",
         "ccr10","cxcr1","cxcr2","cxcr3","cxcr4","cxcr5","cxcr6","cx3cr1",
         "xxcr1","ackr1","ackr2","ackr3","ackr4","ccrl2")
for(i in seq){
  for(j in seq){
    a <- read_csv("data/sequence/gpcr/alignment_csv/ALL_classa_df.csv") %>%
      filter(protein %in% c(i) ) %>%
      dplyr::select(-class, -seq, -protein)
    a <- a %>%
      pivot_longer(cols = 1:ncol(a), names_to = "target_gnccn", values_to = "resid_a")
    b <- read_csv("data/sequence/gpcr/alignment_csv/ALL_classa_df.csv") %>%
      filter(protein %in% c(j) ) %>%
      dplyr::select(-class, -seq, -protein)
    b <- b %>%
      pivot_longer(cols = 1:ncol(b), names_to = "target_gnccn", values_to = "resid_b")
    ab <- left_join(a,b)
    rm(a,b)
    ab <- ab %>%
      mutate(identity = case_when(
        resid_a == resid_b ~ 1,
        resid_a != resid_b ~ 0
      ))
    ab <- ab %>% separate(col = target_gnccn, into = c("a", "target_gnccn"), sep = "gn")
    ab <- ab %>% dplyr::select(-a)
    ab$pair <- c(paste0(i, "_", j))
    # ab <- ab %>% dplyr::select(target_gnccn, identity, pair)
    temp <- rbind(temp, ab)
    rm(ab)
  }
}

temp$pair <- toupper(temp$pair)
write_csv(temp, "data/sequence/gpcr/processed/CKR_ALL_BY_ALL_PAIRWISE_IDENTITY_BY_POS.csv" )
