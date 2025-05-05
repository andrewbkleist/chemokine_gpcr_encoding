# Name:     04_vmip_ackr_comparisons.R
# Updated:  20201122
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggridges)

##### 1: POSITION BASED SCORE COMPARISONS INTERFACE - POINT PLOT ###############

  # (1.1) FUNCTION -----------------------------------------------------------
  GetScorePlotRinCK <- function(POSITION){
    rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") 
    rin <- rin %>% filter(source_gnccn == POSITION) %>%
      select(ck, source_gnccn, cc_cxc_lr_score_ck, cc_cxc_lr_score_sd_ck) %>% unique() %>%
      filter(ck != "cx3cl1")
    
    order <- c("ccl5", "ccl20", "vmipii", "cxcl12", "cxcl8")
    rin$ck <- factor(rin$ck, levels = rev(order))
    
    rin <- rin %>%
      ggplot(aes(cc_cxc_lr_score_ck, ck )) +
      geom_errorbar(aes(xmin=cc_cxc_lr_score_ck-cc_cxc_lr_score_sd_ck, xmax=cc_cxc_lr_score_ck+cc_cxc_lr_score_sd_ck), width=.2) +
      geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
      xlim(0,1) +
      theme_minimal()
    print(rin)
  }
  
  # (1.2) MAKE PLOTS  --------------------------------------------------------
  GetScorePlotRinCK("cxb1.1")
  GetScorePlotRinCK("NTc.Cm1")
  GetScorePlotRinCK("NTc.Cm4")


##### 2: POSITION BASED SCORE COMPARISONS INTERFACE - RIDGE PLOT ###############
  
  # plot all T2-T2 interactions
  read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    #filter(cc_cxc_lr_ck > 0.75 & cc_cxc_lr_ckr > 0.75) %>%
    filter(file %in% pdbs) %>%
    select(file, cc_cxc_lr_ck, cc_cxc_lr_ckr) %>%
    ggplot(aes(cc_cxc_lr_ckr, cc_cxc_lr_ck)) +
    geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
    xlim(0.49,1) +
    ylim(0.49,1) +
    theme_minimal()
  
  # plot spectrum of T2-T2 scores from complexes
  pdbs <- c("6lfo", "6wwz", "ngo", "zheng", "4rws")
  rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
    filter(cc_cxc_lr_ck > 0.75 & cc_cxc_lr_ckr > 0.75) %>%
    filter(file %in% pdbs) %>%
    select(ck, source_gnccn, cc_cxc_lr_score_ck) %>% 
    unique()
  
  order <- c("ccl5", "ccl20", "vmipii", "cxcl12", "cxcl8")
  rin$ck <- factor(rin$ck, levels = rev(order))
  
  rin %>%
    ggplot(aes(x = cc_cxc_lr_score_ck, y = ck)) +
    geom_density_ridges(fill = "steelblue3") +
    theme_minimal() 

##### 3: POSITION BASED SCORE COMPARISONS ALL POSITIONS - POINT PLOT ###########
  
  # (3.1) Get CK (i) T2 positions (ii) @ interface -----------------------------
  
    # import T2
    ck.t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      filter(mean > 0.75)
    
    # import interface
    ck.rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
      select(source_gnccn) %>% unique()
    ck.rin <- ck.rin$source_gnccn
    
    # subset T2 for those at interface...
    ck.t2 <- ck.t2 %>% filter(motif %in% ck.rin)
    # ck.t2 <- ck.t2 %>% filter(motif != "NTc.Cm10" & motif != "NTc.Cm9" & motif != "NTc.Cm8" & motif != "NTc.Cm7")
    ck.t2 <- ck.t2$motif
    rm(ck.rin)
  
  # (3.2) Subset scores by T2 + interface + subset chemokines ----------------
  
    # import scores
    ck.score <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
    ck.vmip <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
    
    ck.score <- bind_rows(ck.score, ck.vmip)
    rm(ck.vmip)    
    
    # subset scores by interface and T2
    ck.score <- ck.score %>% filter(position %in% ck.t2)
    rm(ck.t2)
    
    # choose chemokines of interest
    ck.score <- ck.score %>% filter(protein %in% c("ccl5", "vmip2xhhv8p")) %>%
      select(protein, position, mean, sd) %>% unique()
  
  # (3.3) Order and plot -----------------------------------------------------
    
    # order
    order.ccn <- as.factor(c("NTc.Cm10","NTc.Cm9","NTc.Cm8","NTc.Cm7","NTc.Cm6","NTc.Cm5","NTc.Cm4","NTc.Cm3","NTc.Cm2","NTc.Cm1", "CX.1","CX.2","CX.3","CX.4","CX.5","cxb1.1","cxb1.2","cxb1.3","cxb1.4","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.9","cxb1.10","cxb1.11","cxb1.12","cxb1.13","cxb1.14","cxb1.15","cxb1.16","cxb1.17","cxb1.18","cxb1.19","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.1","b1b2.2","b1b2.3","b1b2.4","b1b2.5","b1b2.6","b1b2.7","b1b2.8","b1b2.9","b1b2.10","b1b2.11","b1b2.12","b1b2.13","b1b2.14","b1b2.15","b1b2.16","b1b2.17","b1b2.18","b1b2.19","b1b2.20","b1b2.21","b1b2.22","b1b2.23","b1b2.24","b1b2.25","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.1","b2b3.2","b2b3.3","b2b3.4","b2b3.5","b2b3.6","b2b3.7","b2b3.8","b2b3.9","b2b3.10","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","b3h.5","b3h.6","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10"))
    ck.score$position <- factor(ck.score$position, levels = order.ccn)
    
    # plot
    ck.score %>%
      ggplot(aes(position, mean, fill = protein)) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
      geom_point(shape = 21,   size = 4, stroke = 0.5) +
      theme_minimal() +
      ylim(0,1) +
      scale_fill_manual(values=c("steelblue4", "red4")) +
      theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0))
    
    # pie chart - add label
    ck.score <- ck.score %>% mutate(cc_cxc = case_when(
      mean > 0.5 ~ "cxc",
      mean < 0.5 ~ "cc"
    ))
    ck.score.sum <- ck.score %>% count(protein, cc_cxc)
    ck.score.sum <- ck.score.sum %>% mutate(pct_cc = n / 34)
    
    # pie chart plot
    ck.score.sum %>%
      ggplot(aes(x = "", y= pct_cc)) +
      geom_bar(width = 1,size = 1, stat="identity", color = "white") +
      coord_polar("y") +
      theme_classic() +
      theme(axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank()) +
      facet_grid(. ~ protein)


##### 4: POSITION BASED SCORE COMPARISONS ALL POSITIONS - POINT PLOT - DARC ####

  # (4.1) Get CKR (i) T2 positions (ii) @ interface ----------------------------
  
    # import T2
    ckr.t2 <- read_csv("03_ckr_seq/output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
      filter(mean > 0.75)
    
    # import interface
    ckr.rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
      select(target_gnccn) %>% unique()
    ckr.rin$gn <- c("gn")
    ckr.rin <- ckr.rin %>% unite(col = target_gnccn, c(gn, target_gnccn), sep = "", remove = TRUE)
    ckr.rin <- ckr.rin$target_gnccn
    
    # subset T2 for those at interface...
    ckr.t2 <- ckr.t2 %>% filter(motif %in% ckr.rin)
    ckr.t2 <- ckr.t2 %>% filter(motif != "gnNTr.Cm26" & motif != "gnNTr.Cm11" & motif != "gnNTr.Cm10" & motif != "gnNTr.Cm6" & motif != "gnNTr.CmCp6")
    ckr.t2 <- ckr.t2$motif
    rm(ckr.rin)

  # (4.2) Subset scores by T2 + interface + subset chemokines ----------------
  
    # import scores
    ckr.score <- read_csv("03_ckr_seq/output/CKR_CLASSA_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
    
    # subset scores by interface and T2
    ckr.score <- ckr.score %>% filter(position %in% ckr.t2)
    rm(ckr.t2)
    
    # choose chemokines of interest
    ckr.score <- ckr.score %>% filter(protein %in% c("ccr5", "ackr1", "cxcr2")) %>%
      select(protein, position, mean, sd) %>% unique()

  # (4.3) Order and plot -----------------------------------------------------
  
    # order
    order.gn <- as.factor(c("gn1x23","gn1x24","gn1x28","gn1x31","gn2x63","gn2x66","gn3x29","gn3x33","gn4x65","gn45x51","gn45x52","gnECL2.Cp6","gn5x36","gn5x40","gn5x43","gn6x58","gn6x62","gn6x64","gn6x65","gn7x27","gn7x34","gn7x35"))
    ckr.score$position <- factor(ckr.score$position, levels = order.gn)
    
    # order by number
    # ck.score <- ck.score %>% filter(protein %in% c("vmip2xhhv8p"))
    # ck.score$position <- factor(ck.score$position, levels = ck.score$position[order(ck.score$mean)])
    
    # plot
    ckr.score %>%
      ggplot(aes(position, mean, fill = protein)) +
      geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2) +
      geom_point(shape = 21,   size = 4, stroke = 0.5) +
      theme_minimal() +
      ylim(0,1) +
      scale_fill_manual(values=c("grey40","steelblue4", "navyblue")) +
      theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0))
    
    # pie chart - add label
    ckr.score <- ckr.score %>% mutate(cc_cxc = case_when(
      mean > 0.5 ~ "cxc",
      mean < 0.5 ~ "cc"
    ))
    ckr.score.sum <- ckr.score %>% count(protein, cc_cxc)
    ckr.score.sum <- ckr.score.sum %>% mutate(pct_cc = n / 22)
  
  
