# Name:     core_vs_nterm_struct_sim.R
# Updated:  20201202
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggrepel)

##### 1: CORE VS. NTERM STRUCT. SIM - CHEMOKINE  ###############################

  # (1.1) Import, label --------------------------------------------------------

  # import rin
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    select(file, source_gnccn, target_gnccn, dom1, dom2) #%>%
    # filter(file != "ngo" & file != "zheng")
  
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
  
  
  # (1.2) Core distance --------------------------------------------------------

  # select chemokine core, spread
  core <- rin %>% filter(dom1 != "NTc") %>% 
    select(file, source_gnccn, target_gnccn) %>% 
    unite(col = rin, 2:3, sep = "_", remove = TRUE)
  core$temp <- c(1)
  core <- core %>% spread(rin, temp, fill = 0)
  
  # matrix and hclust
  core <- as.matrix(core)
  rownames(core) <- core[,1]
  core <- core[,-1]
  
  core <- dist(core, method = "euclidean")
  core <- as.data.frame(as.matrix(core))
  core[lower.tri(core)] <- NA
  core$file1 <- rownames(core)
  core <- core %>% gather(file2, dist, 1:(ncol(core)-1))
  core <- core %>% filter(!is.na(dist))
  
  # (1.2) Nterm distance -------------------------------------------------------
  
  # select nterm, spread
  nterm <- rin %>% filter(dom1== "NTc") %>% 
    select(file, source_gnccn, target_gnccn) %>% 
    unite(col = rin, 2:3, sep = "_", remove = TRUE)
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
  
  # (1.3) Combine, normalize ---------------------------------------------------
  
  # combine
  colnames(core)[3] <- c("dist_core")
  colnames(nterm)[3] <- c("dist_nterm")
  
  master <- left_join(core, nterm)
  master <- master %>% unite(col = file, 1:2, sep = "_", remove = TRUE)
  master <- master %>% filter(dist_core != 0 & dist_nterm != 0)
  
  # normalize
  master <- master %>% mutate(dist_core = dist_core/sqrt(192) ) %>%
    mutate(dist_nterm = dist_nterm/sqrt(131))
  
  # gather
  master <- master %>%
    gather(domain, distance, 2:3)
  
  # (1.4) Stats, plot ----------------------------------------------------------
  
  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  sum <- wilcox.test(distance ~ domain, data = master, exact = FALSE)
  sum$p.value
  
  order <- c("dist_core", "dist_nterm")
  master$domain <- factor(master$domain, levels = rev(order))
  
  master %>%
    ggplot(aes(domain, distance)) +
    geom_violin(trim = FALSE) +
    #geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 0.5) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, colour = "black", fill = "white", stroke = 0.5) +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    #geom_boxplot() +
    theme_minimal() +
    ylim(0,1.2) +
    coord_flip()
  
  
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
  
  
  # (2.2) Core distance --------------------------------------------------------
  
  # select chemokine core, spread
  core <- rin %>% filter(dom2 != "NTr") %>% 
    select(file, source_gnccn, target_gnccn) %>% 
    unite(col = rin, 2:3, sep = "_", remove = TRUE)
  core$temp <- c(1)
  core <- core %>% spread(rin, temp, fill = 0)
  
  # matrix and hclust
  core <- as.matrix(core)
  rownames(core) <- core[,1]
  core <- core[,-1]
  
  core <- dist(core, method = "euclidean")
  core <- as.data.frame(as.matrix(core))
  core[lower.tri(core)] <- NA
  core$file1 <- rownames(core)
  core <- core %>% gather(file2, dist, 1:(ncol(core)-1))
  core <- core %>% filter(!is.na(dist))
  
  # (2.3) Nterm distance -------------------------------------------------------
  
  # select nterm, spread
  nterm <- rin %>% filter(dom2== "NTr") %>% 
    select(file, source_gnccn, target_gnccn) %>% 
    unite(col = rin, 2:3, sep = "_", remove = TRUE)
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
  
  
  # (2.4) ECL2 distance -------------------------------------------------------
  
  # select nterm, spread
  ecl2 <- rin %>% filter(dom2== "ECL2") %>% 
    select(file, source_gnccn, target_gnccn) %>% 
    unite(col = rin, 2:3, sep = "_", remove = TRUE)
  ecl2$temp <- c(1)
  ecl2 <- ecl2 %>% spread(rin, temp, fill = 0)
  
  # matrix and hclust
  ecl2 <- as.matrix(ecl2)
  rownames(ecl2) <- ecl2[,1]
  ecl2 <- ecl2[,-1]
  
  ecl2 <- dist(ecl2, method = "euclidean")
  ecl2 <- as.data.frame(as.matrix(ecl2))
  ecl2[lower.tri(ecl2)] <- NA
  ecl2$file1 <- rownames(ecl2)
  ecl2 <- ecl2 %>% gather(file2, dist, 1:(ncol(ecl2)-1))
  ecl2 <- ecl2 %>% filter(!is.na(dist))
  
  # (2.4) Combine, normalize ---------------------------------------------------
  
  # combine
  colnames(core)[3] <- c("dist_core")
  colnames(nterm)[3] <- c("dist_nterm")
  colnames(ecl2)[3] <- c("dist_ecl2")
  
  master <- left_join(core, nterm)
  master <- left_join(master, ecl2)
  
  master <- master %>% unite(col = file, 1:2, sep = "_", remove = TRUE)
  master <- master %>% filter(dist_core != 0 & dist_nterm != 0 & dist_ecl2 != 0)
  
  # normalize
  master <- master %>% mutate(dist_core = dist_core/sqrt(228)) %>%
    mutate(dist_nterm = dist_nterm/sqrt(95)) %>%
    mutate(dist_ecl2 = dist_ecl2/sqrt(71))
  
  # gather
  master <- master %>%
    gather(domain, distance, 2:4)
  
  # (2.5) Stats, plot ----------------------------------------------------------
  
  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  core.nt <- master %>% filter(domain != "dist_ecl2")
  core.nt <- wilcox.test(distance ~ domain, data = core.nt, exact = FALSE)
  core.nt$p.value
  
  core.ecl2 <- master %>% filter(domain != "dist_nterm")
  core.ecl2 <- wilcox.test(distance ~ domain, data = core.ecl2, exact = FALSE)
  core.ecl2$p.value
  
  nt.ecl2 <- master %>% filter(domain != "dist_core")
  nt.ecl2 <- wilcox.test(distance ~ domain, data = nt.ecl2, exact = FALSE)
  nt.ecl2$p.value
  
  
  order <- c("dist_core", "dist_ecl2", "dist_nterm")
  master$domain <- factor(master$domain, levels = rev(order))
  
  master %>%
    ggplot(aes(domain, distance)) +
    geom_violin(trim = FALSE) +
    #geom_point(shape = 21, colour = "black", fill = "white", size = 2, stroke = 0.5) +
    geom_dotplot(binaxis='y', stackdir='center', dotsize=0.8, colour = "black", fill = "white", stroke = 0.5) +
    stat_summary(fun.y=median, geom="point", shape=23, size=2) +
    #geom_boxplot() +
    theme_minimal() +
    ylim(0,1.2) +
    coord_flip()

    