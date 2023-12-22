source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

temp <- as.data.frame(NULL)
seq <- c("ccr1","ccr2","ccr3","ccr4","ccr5","ccr6","ccr7","ccr8","ccr9",
         "ccr10","cxcr1","cxcr2","cxcr3","cxcr4","cxcr5","cxcr6","cx3cr1",
         "xxcr1","ackr1","ackr2","ackr3","ackr4","ccrl2")
for(i in seq){
  for(j in seq){
    a <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
      filter(protein %in% c(i) ) %>%
      dplyr::select(motif)
    b <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
      filter(protein %in% c(j) ) %>%
      dplyr::select(motif)
    a.only <- dplyr::setdiff(a,b)
    a.only$type <- c("a_not_b")
    b.only <- dplyr::setdiff(b,a)
    b.only$type <- c("b_not_a")
    ab <- dplyr::intersect(a,b)
    ab$type <- c("a_and_b")
    data <- rbind(a.only, ab, b.only)
    data$pair <- c(paste0(i, "_", j))
    rm(a,b,a.only,ab,b.only)
    temp <- rbind(temp, data)
    rm(data)
  }
}

temp$pair <- toupper(temp$pair)
write_csv(temp, "data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CKR_ECL2.csv" )
