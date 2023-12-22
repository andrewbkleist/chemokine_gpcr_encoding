source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
# (1) GENERATE MERS - ECL2A ----------------------------------------------------
  # import alignment
  aln <- readAAMultipleAlignment("data/sequence/gpcr/alignments/ECL2A_CYSLESS_NOGAP.fasta")                   
  aln.df <- as.data.frame(as.matrix(aln))
  
  # 2 MERS ---------------------------------------------------------------------
  # loop over sequences
  mer2 <- NULL
  for (i in 1:nrow(aln.df)){  
    x <- 2
    z <- 1
    for (j in 1:(ncol(aln.df)-1)){
      motif <- aln.df[i, j:x]                        
      colnames(motif) <- c("a", "b")
      motif$file <- rownames(aln.df[i,])            
      motif$motif_no <- z                           
      mer2 <- rbind(mer2, motif)
      z <- z + 1
      x <- x + 1
    }
  }
  rm(i, j, x, z, motif)
  
  # write raw output
  write.csv(mer2, "data/motif/raw/CKR_UNSTRUCTURED_ECL2A_2MERS_RAW_CYSLESS.csv")
  rm(mer2) 
  
  # 3 MERS ---------------------------------------------------------------------
  # loop over sequences
  mer3 <- NULL
  for (i in 1:nrow(aln.df)){  
    x <- 3
    z <- 1
    for (j in 1:(ncol(aln.df)-2)){
      motif <- aln.df[i, j:x]                        
      colnames(motif) <- c("a", "b", "c")
      motif$file <- rownames(aln.df[i,])            
      motif$motif_no <- z                           
      mer3 <- rbind(mer3, motif)
      z <- z + 1
      x <- x + 1
    }
  }
  rm(i, j, x, z, motif)
  
  # write raw output
  write.csv(mer3, "data/motif/raw/CKR_UNSTRUCTURED_ECL2A_3MERS_RAW_CYSLESS.csv")
  rm(mer3)
  
  # 4 MERS ---------------------------------------------------------------------
  # loop over sequences
  mer4 <- NULL
  for (i in 1:nrow(aln.df)){  
    x <- 4
    z <- 1
    for (j in 1:(ncol(aln.df)-3)){
      motif <- aln.df[i, j:x]                        
      colnames(motif) <- c("a", "b", "c", "d")
      motif$file <- rownames(aln.df[i,])            
      motif$motif_no <- z                           
      mer4 <- rbind(mer4, motif)
      z <- z + 1
      x <- x + 1
    }
  }
  rm(i, j, x, z, motif)
  
  # write raw output
  write.csv(mer4, "data/motif/raw/CKR_UNSTRUCTURED_ECL2A_4MERS_RAW_CYSLESS.csv")
  rm(mer4)

    
# (2) GENERATE MERS - ECL2B ----------------------------------------------------
  # import alignment
  aln <- readAAMultipleAlignment("data/sequence/gpcr/alignments/ECL2B_CYSLESS_NOGAP.fasta")                   
  aln.df <- as.data.frame(as.matrix(aln))
  
  # 2 MERS ---------------------------------------------------------------------
  # loop over sequences
  mer2 <- NULL
  for (i in 1:nrow(aln.df)){  
    x <- 2
    z <- 1
    for (j in 1:(ncol(aln.df)-1)){
      motif <- aln.df[i, j:x]                        
      colnames(motif) <- c("a", "b")
      motif$file <- rownames(aln.df[i,])            
      motif$motif_no <- z                           
      mer2 <- rbind(mer2, motif)
      z <- z + 1
      x <- x + 1
    }
  }
  rm(i, j, x, z, motif)
  
  # write raw output
  write.csv(mer2, "data/motif/raw/CKR_UNSTRUCTURED_ECL2B_2MERS_RAW_CYSLESS.csv")
  rm(mer2) 
  
  # 3 MERS ---------------------------------------------------------------------
  # loop over sequences
  mer3 <- NULL
  for (i in 1:nrow(aln.df)){  
    x <- 3
    z <- 1
    for (j in 1:(ncol(aln.df)-2)){
      motif <- aln.df[i, j:x]                        
      colnames(motif) <- c("a", "b", "c")
      motif$file <- rownames(aln.df[i,])            
      motif$motif_no <- z                           
      mer3 <- rbind(mer3, motif)
      z <- z + 1
      x <- x + 1
    }
  }
  rm(i, j, x, z, motif)
  
  # write raw output
  write.csv(mer3, "data/motif/raw/CKR_UNSTRUCTURED_ECL2B_3MERS_RAW_CYSLESS.csv")
  rm(mer3)
  
  # 4 MERS ---------------------------------------------------------------------
  # loop over sequences
  mer4 <- NULL
  for (i in 1:nrow(aln.df)){  
    x <- 4
    z <- 1
    for (j in 1:(ncol(aln.df)-3)){
      motif <- aln.df[i, j:x]                        
      colnames(motif) <- c("a", "b", "c", "d")
      motif$file <- rownames(aln.df[i,])            
      motif$motif_no <- z                           
      mer4 <- rbind(mer4, motif)
      z <- z + 1
      x <- x + 1
    }
  }
  rm(i, j, x, z, motif)
  
  # write raw output
  write.csv(mer4, "data/motif/raw/CKR_UNSTRUCTURED_ECL2B_4MERS_RAW_CYSLESS.csv")
  rm(mer4)
  
# (3) TIDY, ADD INFORMATION, & WRITE OUTPUT - ECL2A ----------------------------
  mer2 <- Raw2Motif("data/motif/raw/CKR_UNSTRUCTURED_ECL2A_2MERS_RAW_CYSLESS.csv",
                    "data/lookup/cc_cxc_ack.csv", "mer2", "none")
  
  mer3 <- Raw2Motif("data/motif/raw/CKR_UNSTRUCTURED_ECL2A_3MERS_RAW_CYSLESS.csv",
                    "data/lookup/cc_cxc_ack.csv", "mer3", "none")
  
  mer4 <- Raw2Motif("data/motif/raw/CKR_UNSTRUCTURED_ECL2A_4MERS_RAW_CYSLESS.csv",
                    "data/lookup/cc_cxc_ack.csv", "mer4", "none")
  
  
  # for 3- and 4- mers, add masks
  mer3.mask1 <- mer3
  mer3.mask1$b <- c("x")
  mer3.mask1$mask <- c("B")
  mer3 <- bind_rows(mer3, mer3.mask1)
  rm(mer3.mask1)
  
  mer4.mask1 <- mer4
  mer4.mask1$b <- c("x")
  mer4.mask1$mask <- c("B")
  
  mer4.mask2 <- mer4
  mer4.mask2$c <- c("x")
  mer4.mask2$mask <- c("C")
  
  mer4.mask3 <- mer4
  mer4.mask3$b <- c("x")
  mer4.mask3$c <- c("x")
  mer4.mask3$mask <- c("BC")
  mer4 <- bind_rows(mer4, mer4.mask1, mer4.mask2, mer4.mask3)
  rm(mer4.mask1, mer4.mask2, mer4.mask3)
  
  # make into "words"
  mer2 <- unite(mer2, col = motif, 6:7,  sep = "")
  mer3 <- unite(mer3, col = motif, 6:8,  sep = "")
  mer4 <- unite(mer4, col = motif, 6:9,  sep = "")
  
  # bind
  data <- bind_rows(mer2, mer3, mer4)
  
  # add unique seqid for each 951 sequences
  seqid <- read_csv("data/motif/raw/CKR_UNSTRUCTURED_ECL2A_2MERS_RAW_CYSLESS.csv")
  seqid <- seqid %>% select(file) %>% distinct()
  seqid$seqid <- c(1:nrow(seqid))
  data$seqid <- seqid$seqid[match(unlist(data$file), seqid$file)]
  
  # write edited version of 3 mer
  # write_csv(data, "data/motif/processed/CKR_UNSTRUCTURED_ECL2A_ALL_MERS_CYSLESS.csv")
  # WRITTEN 2023.11.27
  rm(data, mer2, mer3, mer4)

    
# (3) TIDY, ADD INFORMATION, & WRITE OUTPUT - ECL2B ----------------------------
  mer2 <- Raw2Motif("data/motif/raw/CKR_UNSTRUCTURED_ECL2B_2MERS_RAW_CYSLESS.csv",
                    "data/lookup/cc_cxc_ack.csv", "mer2", "none")
  
  mer3 <- Raw2Motif("data/motif/raw/CKR_UNSTRUCTURED_ECL2B_3MERS_RAW_CYSLESS.csv",
                    "data/lookup/cc_cxc_ack.csv", "mer3", "none")
  
  mer4 <- Raw2Motif("data/motif/raw/CKR_UNSTRUCTURED_ECL2B_4MERS_RAW_CYSLESS.csv",
                    "data/lookup/cc_cxc_ack.csv", "mer4", "none")
  
  
  # for 3- and 4- mers, add masks
  mer3.mask1 <- mer3
  mer3.mask1$b <- c("x")
  mer3.mask1$mask <- c("B")
  mer3 <- bind_rows(mer3, mer3.mask1)
  rm(mer3.mask1)
  
  mer4.mask1 <- mer4
  mer4.mask1$b <- c("x")
  mer4.mask1$mask <- c("B")
  
  mer4.mask2 <- mer4
  mer4.mask2$c <- c("x")
  mer4.mask2$mask <- c("C")
  
  mer4.mask3 <- mer4
  mer4.mask3$b <- c("x")
  mer4.mask3$c <- c("x")
  mer4.mask3$mask <- c("BC")
  mer4 <- bind_rows(mer4, mer4.mask1, mer4.mask2, mer4.mask3)
  rm(mer4.mask1, mer4.mask2, mer4.mask3)
  
  # make into "words"
  mer2 <- unite(mer2, col = motif, 6:7,  sep = "")
  mer3 <- unite(mer3, col = motif, 6:8,  sep = "")
  mer4 <- unite(mer4, col = motif, 6:9,  sep = "")
  
  # bind
  data <- bind_rows(mer2, mer3, mer4)
  
  # add unique seqid for each 951 sequences
  seqid <- read_csv("data/motif/raw/CKR_UNSTRUCTURED_ECL2B_2MERS_RAW_CYSLESS.csv")
  seqid <- seqid %>% select(file) %>% distinct()
  seqid$seqid <- c(1:nrow(seqid))
  data$seqid <- seqid$seqid[match(unlist(data$file), seqid$file)]
  
  
  # write edited version of 3 mer
  # write_csv(data, "data/motif/processed/CKR_UNSTRUCTURED_ECL2B_ALL_MERS_CYSLESS.csv")
  # WRITTEN 2023.11.27
  rm(data, mer2, mer3, mer4)
  