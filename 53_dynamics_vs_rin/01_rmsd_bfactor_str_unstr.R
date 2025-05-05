# Name:     02_rmsd_bfactor_str_unstr.R
# Updated:  20210512
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(bio3d)

##### 1: CHEMOKINE - NMR ENSEMBLE RMSD #########################################

  # GetRMSF <- function(INPUT.FILE, LOOKUP.FILE, )
  # read pdb
  pdb <- read.pdb("51_dynamics_vs_rin/data/raw/2kec.pdb", multi=TRUE)
  
  # select, fit CA
  ca.inds <- atom.select(pdb, elety="CA")
  xyz <- fit.xyz(fixed=pdb$xyz[1,], mobile=pdb$xyz[2:nrow(pdb$xyz),], # recall that pdb$xyz 
                 fixed.inds=ca.inds$xyz,                   # rows are ensemble nos (1-20)
                 mobile.inds=ca.inds$xyz)                  # and cols are atom 1 xyz, atom 2 xyz...
  
      # rd <- rmsd(xyz[1,ca.inds$xyz], xyz[,ca.inds$xyz]) # RMSD each frame to first (?), entire protein
      # rd <- rmsd(pdb$xyz, fit=TRUE) # full structure RMSD of NMR ensemble
      # plot(rd, typ="l", ylab="RMSD", xlab="Frame No.")
  
  # RMSF calculation
  rf <- rmsf(xyz[,ca.inds$xyz])
      # plot(rf, ylab="RMSF", xlab="Residue Position", typ="l")
  rf <- as.data.frame(rf)
  rf$resno <- rownames(rf)
  # rm(xyz, ca.inds, pdb)
  
  # map generic numbering
  lookup <- read_csv("51_dynamics_vs_rin/data/processed/lookup_resno_to_generic.csv")
  rf$gn <- lookup$gn[match(unlist(rf$resno), lookup$resno)]
  rm(lookup)

  # map tp pdb
  temp <- pdb$atom
  temp$b <- rf$rf[match(unlist(temp$resno), rf$resno)]
  
  write.pdb(pdb=pdb, b = temp$b, file = "51_dynamics_vs_rin/output/cxcl12_rmsf_ca.pdb")
  
  # order, write
  rf <- rf %>% select(resno, gn, rf)
  write_csv(rf, "51_dynamics_vs_rin/data/processed/chemokine_rmsf.csv")
  
  
  
##### 2: CHEMOKINE - MAP CONTACT INFORMATION #########################################
  
  # import
  rf <- read_csv("51_dynamics_vs_rin/data/processed/chemokine_rmsf.csv")
  rf <- rf %>%  mutate(sse = case_when(
    grepl("NT", rf["gn"]) ~ "NT",
    !grepl("NT", rf["gn"]) ~ "core"
  ))
  
  
  
  str <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  
  # subset, join
  str <- str %>% select(source_gnccn, all_para_ck, no_contacts_all_ck) %>% unique()
  colnames(str)[1] <- c("gn")
  rf <- left_join(rf, str)
  rf <- rf %>% filter(!is.na(no_contacts_all_ck))
  rm(str)
  
  # plot
  rf %>%
    # filter(rf < 4) %>%
    ggplot(aes(all_para_ck, no_contacts_all_ck)) +
    geom_jitter() +
    # xlim(0,9) +
    # ylim(0,25) +
    # geom_smooth() +
    theme_minimal() 
  
  # p + geom_bin2d(binwidth=c(1, 2))
  
  
  # import
  rf <- read_csv("50_dynamics_vs_rin/data/processed/chemokine_rmsf.csv")
  str <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
  # seq <- read_csv("05_integrate/output/CK_CONS_CCCXC_SNP_CAN.csv")

  
  # 
  str <- str %>% select(source_gnccn, all_para_ck, no_contacts_all_ck) %>% unique()
  colnames(str)[1] <- c("gn")
  
  test <- left_join(rf, str)
  test <- test %>% filter(!is.na(no_contacts_all_ck))
  
  model.wt <- SummarizeGrowth(test$rf, test$no_contacts_all_ck)
  df.predicted <- data.frame(rf = test$rf, pred.wt = predict(model.wt$model))
  
  pl <- test %>%
    # filter(rf < 4) %>%
    ggplot(aes(rf, no_contacts_all_ck)) +
    geom_jitter() +
     # xlim(0,9) +
     # ylim(0,25) +
    theme_minimal() 
  
  
  df.predicted <- data.frame(rf = test$rf, pred.wt = predict(model.wt$model))
  pl + geom_line(data=df.predicted, aes(y=df.predicted$pred.wt), color="red")
  
  
  cor.test(test$no_contacts_all_ck, test$rf, method = "pearson", conf.level = 0.95)
  
  
  
  
  
  
  
  
  
