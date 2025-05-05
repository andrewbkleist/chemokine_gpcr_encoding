# Name:     01_motif_gen.R
# Updated:  20201218
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(Biostrings)
library(seqRFLP)
library(stringr)

# **THIS SCRIPT WAS RUN IN PRIOR TO THE ORIGINAL SUBMISSION (c. 11/2020; 
# see /kleist_2020 foler) AND NOT RE-RUN DURING REVISIONS DUE TO LONG RUN-TIMES
# AND NO CHANGES TO CODE; INPUT/OUTPUT PATHS AND REQUIRED INPUTS WERE MODIFIED
# TO ACCOMODATE NEW FOLDER HIERARCHIES

##### 1: GENERATE MERS - NTERM #################################################
  
  # import alignment
  aln <- readAAMultipleAlignment("08_ckr_motif/data/raw/NTERM_CYSLESS_NOGAP.fasta")                   
  aln.df <- as.data.frame(as.matrix(aln))
  
  
  # 2 MERS ------------------------------------------------------------------------
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
  write.csv(mer2, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_2MERS_RAW_CYSLESS.csv")
  rm(mer2) 
  
  # 3 MERS ------------------------------------------------------------------------
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
  write.csv(mer3, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_3MERS_RAW_CYSLESS.csv")
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
  write.csv(mer4, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_4MERS_RAW_CYSLESS.csv")
  rm(mer4)
  
  
##### 2: GENERATE MERS - ECL2A #################################################
  
  # import alignment
  aln <- readAAMultipleAlignment("08_ckr_motif/data/raw/ECL2A_CYSLESS_NOGAP.fasta")                   
  aln.df <- as.data.frame(as.matrix(aln))
  
  # 2 MERS ------------------------------------------------------------------------
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
  write.csv(mer2, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_2MERS_RAW_CYSLESS.csv")
  rm(mer2) 
  
  # 3 MERS ------------------------------------------------------------------------
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
  write.csv(mer3, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_3MERS_RAW_CYSLESS.csv")
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
  write.csv(mer4, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_4MERS_RAW_CYSLESS.csv")
  rm(mer4)

    
##### 2: GENERATE MERS - ECL2B #################################################
  
  # import alignment
  aln <- readAAMultipleAlignment("08_ckr_motif/data/raw/ECL2B_CYSLESS_NOGAP.fasta")                   
  aln.df <- as.data.frame(as.matrix(aln))
  
  # 2 MERS ------------------------------------------------------------------------
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
  write.csv(mer2, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_2MERS_RAW_CYSLESS.csv")
  rm(mer2) 
  
  # 3 MERS ------------------------------------------------------------------------
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
  write.csv(mer3, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_3MERS_RAW_CYSLESS.csv")
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
  write.csv(mer4, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_4MERS_RAW_CYSLESS.csv")
  rm(mer4)
  
##### 3: TIDY MERS - NTERM #####################################################
  
  raw2motif <- function(FILE, LOOKUP_CLASS, MER, MASK){
    
    # import data
    data <- read_csv(FILE)
    data <- data[,-1]
    
    # extract chemokine name
    name <- strsplit(data$file, "/", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    colnames(name) <- ("protein")
    name <- strsplit(name[,1], "_", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    rownames(name) <- 1:nrow(name)
    colnames(name) <- c("protein")
    name <- as.data.frame(name)
    data <- cbind(data, name)
    rm(name)

    # add classification variables (CC, CXC, ACK, CX3L, XC)
    cc.cxc.ack <- read.csv(LOOKUP_CLASS)
    data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
    rm(cc.cxc.ack)
    
    # add mer designation
    data$mer <- c(MER)
    
    # reorder, filter, remove caps
    if(MER == "mer2"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b) %>%
        filter(a != "-" & b != "-" )
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      data <- data %>% select(-a,-b) %>% bind_cols(a,b) %>% 
        select(protein, class, file, mer, motif_no, a, b)
      
    } else if (MER == "mer3"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b, c) %>%
        filter(a != "-" & b != "-" & c != "-") 
      
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      c <- as.data.frame(toupper(data$c))
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      colnames(c) <- c("c")
      data <- data %>% select(-a,-b,-c) %>% bind_cols(a,b,c) %>% 
        select(protein, class, file, mer, motif_no, a, b, c)
      
    } else if (MER == "mer4"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b, c, d) %>%
        filter(a != "-" & b != "-" & c != "-" & d != "-") 
      
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      c <- as.data.frame(toupper(data$c))
      d <- as.data.frame(toupper(data$d))
      
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      colnames(c) <- c("c")
      colnames(d) <- c("d")
      data <- data %>% select(-a,-b,-c, -d) %>% bind_cols(a,b,c,d) %>% 
        select(protein, class, file, mer, motif_no, a, b, c, d)
    }
    
    data$mask <- c(MASK)
    
    return(data)
    rm(data, a, b, c, d)
  }
  
  mer2 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_4MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer2", "none")
  
  mer3 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_4MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer3", "none")
  
  mer4 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_4MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer4", "none")
  

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
  
  # write edited version of 3 mer
  write_csv(data, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ALL_MERS_CYSLESS.csv")
  rm(data, mer2, mer3, mer4)
  
  
##### 4: TIDY MERS - ECL2A ######################################################
  
  raw2motif <- function(FILE, LOOKUP_CLASS, MER, MASK){
    
    # import data
    data <- read_csv(FILE)
    data <- data[,-1]
    
    # extract chemokine name
    name <- strsplit(data$file, "/", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    colnames(name) <- ("protein")
    name <- strsplit(name[,1], "_", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    rownames(name) <- 1:nrow(name)
    colnames(name) <- c("protein")
    name <- as.data.frame(name)
    data <- cbind(data, name)
    rm(name)
    
    # add classification variables (CC, CXC, ACK, CX3L, XC)
    cc.cxc.ack <- read.csv(LOOKUP_CLASS)
    data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
    rm(cc.cxc.ack)
    
    # add mer designation
    data$mer <- c(MER)
    
    # reorder, filter, remove caps
    if(MER == "mer2"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b) %>%
        filter(a != "-" & b != "-" )
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      data <- data %>% select(-a,-b) %>% bind_cols(a,b) %>% 
        select(protein, class, file, mer, motif_no, a, b)
      
    } else if (MER == "mer3"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b, c) %>%
        filter(a != "-" & b != "-" & c != "-") 
      
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      c <- as.data.frame(toupper(data$c))
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      colnames(c) <- c("c")
      data <- data %>% select(-a,-b,-c) %>% bind_cols(a,b,c) %>% 
        select(protein, class, file, mer, motif_no, a, b, c)
      
    } else if (MER == "mer4"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b, c, d) %>%
        filter(a != "-" & b != "-" & c != "-" & d != "-") 
      
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      c <- as.data.frame(toupper(data$c))
      d <- as.data.frame(toupper(data$d))
      
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      colnames(c) <- c("c")
      colnames(d) <- c("d")
      data <- data %>% select(-a,-b,-c, -d) %>% bind_cols(a,b,c,d) %>% 
        select(protein, class, file, mer, motif_no, a, b, c, d)
    }
    
    data$mask <- c(MASK)
    
    return(data)
    rm(data, a, b, c, d)
  }
  
  mer2 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_2MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer2", "none")
  
  mer3 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_3MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer3", "none")
  
  mer4 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_4MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer4", "none")
  
  
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
  seqid <- read_csv("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_2MERS_RAW_CYSLESS.csv")
  seqid <- seqid %>% select(file) %>% distinct()
  seqid$seqid <- c(1:nrow(seqid))
  data$seqid <- seqid$seqid[match(unlist(data$file), seqid$file)]
  
  # write edited version of 3 mer
  write_csv(data, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2A_ALL_MERS_CYSLESS.csv")
  rm(data, mer2, mer3, mer4)

    
##### 5: TIDY MERS - ECL2B ######################################################
  
  raw2motif <- function(FILE, LOOKUP_CLASS, MER, MASK){
    
    # import data
    data <- read_csv(FILE)
    data <- data[,-1]
    
    # extract chemokine name
    name <- strsplit(data$file, "/", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    colnames(name) <- ("protein")
    name <- strsplit(name[,1], "_", fixed = TRUE)
    name <- t(as.data.frame(lapply(name, "[", 1)))
    rownames(name) <- 1:nrow(name)
    colnames(name) <- c("protein")
    name <- as.data.frame(name)
    data <- cbind(data, name)
    rm(name)
    
    # add classification variables (CC, CXC, ACK, CX3L, XC)
    cc.cxc.ack <- read.csv(LOOKUP_CLASS)
    data$class <- cc.cxc.ack$class[match(unlist(data$protein), cc.cxc.ack$ck)]
    rm(cc.cxc.ack)
    
    # add mer designation
    data$mer <- c(MER)
    
    # reorder, filter, remove caps
    if(MER == "mer2"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b) %>%
        filter(a != "-" & b != "-" )
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      data <- data %>% select(-a,-b) %>% bind_cols(a,b) %>% 
        select(protein, class, file, mer, motif_no, a, b)
      
    } else if (MER == "mer3"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b, c) %>%
        filter(a != "-" & b != "-" & c != "-") 
      
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      c <- as.data.frame(toupper(data$c))
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      colnames(c) <- c("c")
      data <- data %>% select(-a,-b,-c) %>% bind_cols(a,b,c) %>% 
        select(protein, class, file, mer, motif_no, a, b, c)
      
    } else if (MER == "mer4"){
      
      data <- data %>% select(protein, class, file, mer, motif_no, a, b, c, d) %>%
        filter(a != "-" & b != "-" & c != "-" & d != "-") 
      
      # force all caps
      a <- as.data.frame(toupper(data$a))
      b <- as.data.frame(toupper(data$b))
      c <- as.data.frame(toupper(data$c))
      d <- as.data.frame(toupper(data$d))
      
      colnames(a) <- c("a")
      colnames(b) <- c("b")
      colnames(c) <- c("c")
      colnames(d) <- c("d")
      data <- data %>% select(-a,-b,-c, -d) %>% bind_cols(a,b,c,d) %>% 
        select(protein, class, file, mer, motif_no, a, b, c, d)
    }
    
    data$mask <- c(MASK)
    
    return(data)
    rm(data, a, b, c, d)
  }
  
  mer2 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_2MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer2", "none")
  
  mer3 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_3MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer3", "none")
  
  mer4 <- raw2motif("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_4MERS_RAW_CYSLESS.csv",
                    "08_ckr_motif/data/processed/cc_cxc_ack.csv", "mer4", "none")
  
  
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
  seqid <- read_csv("08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_2MERS_RAW_CYSLESS.csv")
  seqid <- seqid %>% select(file) %>% distinct()
  seqid$seqid <- c(1:nrow(seqid))
  data$seqid <- seqid$seqid[match(unlist(data$file), seqid$file)]
  
  
  # write edited version of 3 mer
  write_csv(data, "08_ckr_motif/data/processed/CKR_UNSTRUCTURED_ECL2B_ALL_MERS_CYSLESS.csv")
  rm(data, mer2, mer3, mer4)
  