# Name:     03_nterm_struct_sim_length.R
# Updated:  20201207
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggrepel)

##### 1: CORE VS. NTERM STRUCT. SIM - CHEMOKINE  ###############################

  # (1.1) Import, label --------------------------------------------------------

  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    select(file, source_gnccn, target_gnccn, dom1, dom2)# %>%
    # filter(file == "4xt1" | file == "5wb2") %>%
   #filter(file != "ngo" & file != "zheng")
  
  # add labels
  rin <- rin %>% 
    mutate(str_ck = case_when(
      dom1 == "NTc" ~ "NTc",
      dom1 != "NTc" ~ "core_c"
    )) %>% mutate(str_ckr = case_when(
    dom2 == "NTr" ~ "NTr",
    dom2 == "ECL2" ~ "ECL2",
    dom2 != "NTr" & dom2 != "ECL2"  ~ "core_r"
  ))
  
  
  # (1.2) distance function ----------------------------------------------------
  
  
  GetDistance <- function(CCN){
    nterm <- rin %>% filter(source_gnccn == CCN) %>% 
      select(file, source_gnccn, target_gnccn) %>% 
      unite(col = rin, 2:3, sep = "_", remove = TRUE)
    getsize <- nrow(nterm)
    nterm$temp <- c(1)
    nterm <- nterm %>% spread(rin, temp, fill = 0)
    
    # matrix and hclust
    nterm <- as.matrix(nterm)
    rownames(nterm) <- nterm[,1]
    nterm <- nterm[,-1]
    
    nterm <- dist(nterm, method = "euclidean")
    nterm <- as.data.frame(as.matrix(nterm))
    nterm[lower.tri(nterm)] <- NA
    nterm$file1 <- rownames(nterm)
    nterm <- nterm %>% gather(file2, dist, 1:(ncol(nterm)-1))
    nterm <- nterm %>% filter(!is.na(dist))
    
    # normalize
    nterm <- nterm %>% mutate(dist = dist/sqrt(getsize) )
    
    # add label
    nterm$class <- c(CCN)
    return(nterm)
  }
  
  # (1.3) Get distances --------------------------------------------------------
  
  NTc.Cm1 <- GetDistance("NTc.Cm1")
  NTc.Cm2 <- GetDistance("NTc.Cm2")
  NTc.Cm3 <- GetDistance("NTc.Cm3")
  NTc.Cm4 <- GetDistance("NTc.Cm4")
  NTc.Cm5 <- GetDistance("NTc.Cm5")
  NTc.Cm6 <- GetDistance("NTc.Cm6")
  NTc.Cm7 <- GetDistance("NTc.Cm7")
  NTc.Cm8 <- GetDistance("NTc.Cm8")
  NTc.Cm9 <- GetDistance("NTc.Cm9")
  NTc.Cm10 <- GetDistance("NTc.Cm10")
  
  master <- bind_rows(NTc.Cm1, NTc.Cm2, NTc.Cm3, NTc.Cm4, NTc.Cm5,
                      NTc.Cm6, NTc.Cm7, NTc.Cm8, NTc.Cm9, NTc.Cm10)
  
  master <- master %>% filter(file1 != file2)
  
  rm(NTc.Cm1, NTc.Cm2, NTc.Cm3, NTc.Cm4, NTc.Cm5,
     NTc.Cm6, NTc.Cm7, NTc.Cm8, NTc.Cm9, NTc.Cm10)
  
  # (1.4) Plot -----------------------------------------------------------------
  
  order <- c("NTc.Cm1", "NTc.Cm2", "NTc.Cm3", "NTc.Cm4", "NTc.Cm5",
             "NTc.Cm6", "NTc.Cm7", "NTc.Cm8", "NTc.Cm9", "NTc.Cm10")
  master$class <- factor(master$class, levels = rev(order))
  
  master %>%
    #filter(file1 == "4xt1" & file2 == "5wb2") %>%
    #filter(file1 == "5uiw" & file2 == "6wwz") %>%
    #filter(file1 == "4rws" & file2 == "5uiw") %>%
    
    ggplot(aes(class, dist)) +
    geom_violin(trim = FALSE) +
    #geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 0.5) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, colour = "black", fill = "white", stroke = 0.5) +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    #geom_boxplot() +
    theme_minimal() +
    ylim(0,1.2) +
    coord_flip()
  
  #rm(rin, master, order, GetDistance)  
  
##### 2: CORE VS. NTERM STRUCT. SIM - RECEPTOR  ################################
  
  # (2.1) Import, label --------------------------------------------------------
  
  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    select(file, source_gnccn, target_gnccn, dom1, dom2) %>%
    filter(file != "ngo" & file != "zheng")
  
  # add labels
  rin <- rin %>% 
    mutate(str_ck = case_when(
      dom1 == "NTc" ~ "NTc",
      dom1 != "NTc" ~ "core_c"
    )) %>% mutate(str_ckr = case_when(
      dom2 == "NTr" ~ "NTr",
      dom2 == "ECL2" ~ "ECL2",
      dom2 != "NTr" & dom2 != "ECL2"  ~ "core_r"
    ))
  
  
  # (2.2) distance function ----------------------------------------------------
  
  GetDistance <- function(GN){
    nterm <- rin %>% filter(target_gnccn == GN) %>% 
      select(file, source_gnccn, target_gnccn) %>% 
      unite(col = rin, 2:3, sep = "_", remove = TRUE)
    getsize <- nrow(nterm)
    nterm$temp <- c(1)
    nterm <- nterm %>% spread(rin, temp, fill = 0)
    
    # matrix and hclust
    nterm <- as.matrix(nterm)
    rownames(nterm) <- nterm[,1]
    nterm <- nterm[,-1]
    
    nterm <- dist(nterm, method = "euclidean")
    nterm <- as.data.frame(as.matrix(nterm))
    nterm[lower.tri(nterm)] <- NA
    nterm$file1 <- rownames(nterm)
    nterm <- nterm %>% gather(file2, dist, 1:(ncol(nterm)-1))
    nterm <- nterm %>% filter(!is.na(dist))
    
    # normalize
    nterm <- nterm %>% mutate(dist = dist/sqrt(getsize) )
    
    # add label
    nterm$class <- c(GN)
    return(nterm)
  }
  
  # (2.3) Get distances --------------------------------------------------------
  
  NTr.Cm1 <- GetDistance("NTr.Cm1")
  NTr.Cm2 <- GetDistance("NTr.Cm2")
  NTr.Cm3 <- GetDistance("NTr.Cm3")
  NTr.Cm4 <- GetDistance("NTr.Cm4")
  NTr.Cm5 <- GetDistance("NTr.Cm5")
  NTr.Cm6 <- GetDistance("NTr.Cm6")
  NTr.Cm7 <- GetDistance("NTr.Cm7")
  NTr.Cm8 <- GetDistance("NTr.Cm8")
  NTr.Cm9 <- GetDistance("NTr.Cm9")
  NTr.Cm10 <- GetDistance("NTr.Cm10")
  
  master <- bind_rows(NTr.Cm1, NTr.Cm2, NTr.Cm3, NTr.Cm4, NTr.Cm5,
                      NTr.Cm6, NTr.Cm7, NTr.Cm8, NTr.Cm9, NTr.Cm10)
  
  master <- master %>% filter(file1 != file2)
  
  rm(NTr.Cm1, NTr.Cm2, NTr.Cm3, NTr.Cm4, NTr.Cm5,
     NTr.Cm6, NTr.Cm7, NTr.Cm8, NTr.Cm9, NTr.Cm10)
  
  # (2.4) Plot -----------------------------------------------------------------
  
  order <- c("NTr.Cm1", "NTr.Cm2", "NTr.Cm3", "NTr.Cm4", "NTr.Cm5",
             "NTr.Cm6", "NTr.Cm7", "NTr.Cm8", "NTr.Cm9", "NTr.Cm10")
  master$class <- factor(master$class, levels = order)
  
  master %>%
    ggplot(aes(class, dist)) +
    geom_violin(trim = FALSE) +
    #geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 0.5) +
    #geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, colour = "black", fill = "white", stroke = 0.5) +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    #geom_boxplot() +
    theme_minimal() +
    ylim(0,1.2) +
    coord_flip()

  rm(rin, master, order, GetDistance)  