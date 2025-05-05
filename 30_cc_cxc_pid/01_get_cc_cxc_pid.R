# Name:     01_get_cc_cxc_pid.R
# Updated:  20221109
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

##### 1: GET ALL-BY-ALL (AMONG CC/CXC) SEQ IDENTITY ############################

  # (1.1) CHEMOKINE ------------------------------------------------------------
  
  # import, names
  pid.ck <- read.table("30_cc_cxc_pid/data/raw/ck_pid_matrix_mod.txt", header = FALSE)

  names <- pid.ck[,1]
  rownames(pid.ck) <- pid.ck[,1]
  pid.ck <- pid.ck[,-1]
  colnames(pid.ck) <- rownames(pid.ck)
  
  # remove duplicates
  pid.ck <- as.matrix(pid.ck)
  pid.ck[lower.tri(pid.ck, diag = TRUE)] <- 0
  pid.ck <- as.data.frame(pid.ck)

  # gather, unite
  pid.ck$ck1 <- names
  pid.ck <- gather(pid.ck, key = "ck2", value = "pid", 1:46)
  pid.ck <- pid.ck %>% filter(pid != 0)
  
  # cc-cc vs. cxc-cxc vs. cc-cxc, etc
  pid.ck <- pid.ck %>% mutate(type = case_when(
    grepl("CCL", pid.ck$ck1) & grepl("CCL", pid.ck$ck2) ~ "cc_cc",
    grepl("CXCL", pid.ck$ck1) & grepl("CXCL", pid.ck$ck2) ~ "cxc_cxc",
    grepl("CCL", pid.ck$ck1) & grepl("CXCL", pid.ck$ck2) ~ "inter",
    grepl("CXCL", pid.ck$ck1) & grepl("CCL", pid.ck$ck2) ~ "inter"
  ))
  
  # remove CX3CL/XCL
  pid.ck <- pid.ck %>% 
    filter(ck1 !="XCL1") %>%  
    filter(ck1 !="XCL2") %>%  
    filter(ck1 !="CX3CL1") %>%
    filter(ck2 !="XCL1") %>%  
    filter(ck2 !="XCL2") %>%  
    filter(ck2 !="CX3CL1")
  
  # order
  pid.ck$type <- factor(pid.ck$type, levels = c("inter", "cxc_cxc", "cc_cc"))
  
  # plot
  pid.ck %>%
    ggplot(aes(type, pid)) +
    #geom_point() +
    geom_violin(alpha = 0.5, trim = TRUE) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
    ylim(0,100) +
    coord_flip() +
    theme_minimal()
  
  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  stat.ck <- pid.ck %>% filter(type == "inter" | type == "cc_cc")
  stat.ck <- wilcox.test(pid ~ type, data = stat.ck, exact = FALSE)
  stat.ck$p.value
  
  stat.ck <- pid.ck %>% filter(type == "inter" | type == "cxc_cxc")
  stat.ck <- wilcox.test(pid ~ type, data = stat.ck, exact = FALSE)
  stat.ck$p.value
  
  stat.ck <- pid.ck %>% filter(type == "cxc_cxc" | type == "cc_cc")
  stat.ck <- wilcox.test(pid ~ type, data = stat.ck, exact = FALSE)
  stat.ck$p.value
  
  # remove
  rm(pid.ck, stat.ck, names)
  
  # (1.2) RECEPTOR ------------------------------------------------------------
  # import, names
  pid.ckr <- read.table("30_cc_cxc_pid/data/raw/ckr_pid_matrix_mod.csv", header = FALSE)
    
  names <- pid.ckr[,1]
  rownames(pid.ckr) <- pid.ckr[,1]
  pid.ckr <- pid.ckr[,-1]
  colnames(pid.ckr) <- rownames(pid.ckr)
  
  # remove duplicates
  pid.ckr <- as.matrix(pid.ckr)
  pid.ckr[lower.tri(pid.ckr, diag = TRUE)] <- 0
  pid.ckr <- as.data.frame(pid.ckr)
  
  # gather, unite
  pid.ckr$ckr1 <- names
  pid.ckr <- gather(pid.ckr, key = "ckr2", value = "pid", 1:23)
  pid.ckr <- pid.ckr %>% filter(pid != 0)
  
  # remove ACKR etc
  pid.ckr <- pid.ckr 
  pid.ckr <- pid.ckr %>% filter(!grepl("ackr", pid.ckr$ckr1))
  pid.ckr <- pid.ckr %>% filter(!grepl("ackr", pid.ckr$ckr2))
  pid.ckr <- pid.ckr %>% filter(ckr1 != "xcr1_human")
  pid.ckr <- pid.ckr %>% filter(ckr2 != "xcr1_human")
  pid.ckr <- pid.ckr %>% filter(ckr1 != "cx3cr1_human") 
  pid.ckr <- pid.ckr %>% filter(ckr2 != "cx3cr1_human") 
  pid.ckr <- pid.ckr %>% filter(ckr1 != "ccrl2_human") 
  pid.ckr <- pid.ckr %>% filter(ckr2 != "ccrl2_human")
  pid.ckr <- pid.ckr %>% filter(ckr1 != "cx3c1_human") 
  pid.ckr <- pid.ckr %>% filter(ckr2 != "cx3c1_human")
    
  # cc-cc vs. cxc-cxc vs. cc-cxc, etc
  pid.ckr <- pid.ckr %>% mutate(type = case_when(
    grepl("ccr", pid.ckr$ckr1) & grepl("ccr", pid.ckr$ckr2) ~ "cc_cc",
    grepl("cxcr", pid.ckr$ckr1) & grepl("cxcr", pid.ckr$ckr2) ~ "cxc_cxc",
    grepl("ccr", pid.ckr$ckr1) & grepl("cxcr", pid.ckr$ckr2) ~ "inter",
    grepl("cxcr", pid.ckr$ckr1) & grepl("ccr", pid.ckr$ckr2) ~ "inter"
  ))
  
  # order
  pid.ckr$type <- factor(pid.ckr$type, levels = c("inter", "cxc_cxc", "cc_cc"))
  
  # plot
  pid.ckr %>%
    ggplot(aes(type, pid)) +
    #geom_point() +
    geom_violin(alpha = 0.5) +
    stat_summary(fun.y=mean, geom="point", shape=23, size=2) +
    ylim(0,100) +
    coord_flip() +
    theme_minimal()
  
  # summary stats and p-value
  # see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
  stat.ckr <- pid.ckr %>% filter(type == "inter" | type == "cc_cc")
  stat.ckr <- wilcox.test(pid ~ type, data = stat.ckr, exact = FALSE)
  stat.ckr$p.value
  
  stat.ckr <- pid.ckr %>% filter(type == "inter" | type == "cxc_cxc")
  stat.ckr <- wilcox.test(pid ~ type, data = stat.ckr, exact = FALSE)
  stat.ckr$p.value
  
  stat.ckr <- pid.ckr %>% filter(type == "cxc_cxc" | type == "cc_cc")
  stat.ckr <- wilcox.test(pid ~ type, data = stat.ckr, exact = FALSE)
  stat.ckr$p.value
  
  # remove
  rm(pid.ckr, stat.ckr, names)
  
  