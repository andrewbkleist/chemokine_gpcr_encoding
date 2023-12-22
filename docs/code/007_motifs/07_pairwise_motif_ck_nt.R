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
    a <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
      filter(protein %in% c(i) ) %>%
      dplyr::select(motif)
    b <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
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
write_csv(temp, "data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CK_NT.csv" )
