source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

temp <- as.data.frame(NULL)
seq <- c("ccl1","ccl2","ccl3","ccl3l1","ccl4","ccl4l1","ccl5","ccl7","ccl8",
         "ccl11","ccl13","ccl14","ccl15","ccl16","ccl17","ccl18","ccl19",
         "ccl20","ccl21","ccl22","ccl23","ccl24","ccl25","ccl26","ccl27",
         "ccl28","cxcl1","cxcl2","cxcl3","cxcl4","cxcl4l1","cxcl5","cxcl6",
         "cxcl7","cxcl8","cxcl9","cxcl10","cxcl11","cxcl12","cxcl13","cxcl14",
         "cxcl16","cxcl17","cx3cl1","xcl1","xcl2")
for(i in seq){
  for(j in seq){
    a <- read_csv("data/sequence/chemokine/alignment_csv/ALL_para_df.csv") %>%
      filter(protein %in% c(i) ) %>%
      dplyr::select(-class, -seq, -protein)
    a <- a %>%
      pivot_longer(cols = 1:ncol(a), names_to = "source_gnccn", values_to = "resid_a")
    b <- read_csv("data/sequence/chemokine/alignment_csv/ALL_para_df.csv") %>%
      filter(protein %in% c(j) ) %>%
      dplyr::select(-class, -seq, -protein)
    b <- b %>%
      pivot_longer(cols = 1:ncol(b), names_to = "source_gnccn", values_to = "resid_b")
    ab <- left_join(a,b)
    rm(a,b)
    ab <- ab %>%
      mutate(identity = case_when(
        resid_a == resid_b ~ 1,
        resid_a != resid_b ~ 0
      ))
    ab$pair <- c(paste0(i, "_", j))
    # ab <- ab %>% dplyr::select(target_gnccn, identity, pair)
    temp <- rbind(temp, ab)
    rm(ab)
  }
}

temp$pair <- toupper(temp$pair)
write_csv(temp, "data/sequence/chemokine/processed/CK_ALL_BY_ALL_PAIRWISE_IDENTITY_BY_POS.csv" )
