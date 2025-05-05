# Name:     1_prep_sat_mut.R
# Updated:  20210329
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

# NOTE 1 (20210330):
# Downloaded data files manually "saved"cleaned" by removing annotations leaving 
# only raw data, saved as /processed/cxcr4_clean.csv
#
#
# NOTE 2 (20210329):
# Emailed Erik Procko (dataset author) early September; he responded 9/8/2021
# We asked him to clarify "CXCL12 binding conservation score" and "Binding
# conservation score" from Fig 3C. 
#
# He responded that analysis was done as follows:
# "The enrichment ratios for all four expression experiments (anti-myc-FITC and 
# anti-myc-Alexa647) in linear form, i.e. 2^(log2 enrichment ratio), are averaged 
# from all the experiments and only then converted to the log2 ratios. By 
# combining the raw enrichment ratios rather than their log forms, it becomes 
# equivalent to if one had combined all the samples together and treated them 
# as a single pooled sample.  Positional conservation scores were then calculated 
# by averaging the log2 enrichment ratios for all mutations including stop 
# codons.  I'm guessing the reasoning for including stop codons is that the 
# expression tag is at the N-terminus and truncated proteins may still be 
# expressed.  In all of our recent work we ignore the stop codons because 
# they don't give expressed protein.  Wild type is 0."
#
# He also warned to average over all mutations:
# Whenever we average replicate data sets we treat it as combining the replicates 
# to give a single pooled sample, so averaging is with the raw ratios.  We only 
#  average the log2 ratios when calculating a conservation score, which is an 
# abstract term to communicate whether a position is mutationally tolerant or not.
# ...I very much view the CCR5/CXCR4 data as more predictive and qualitative"
#
#
# NOTE 3 (20210330):
# "NA" values reflect missing values: "Missing mutations (<10 reads in the naive
# library) are black." (Fig 1)

##### 1: PREPARE DATA  ###########################################################
  
  # (1.1) IMPORT, ADD RESNO, RESID ---------------------------------------------
  # import data
  data <- read_csv("09_sat_mut_data/data/processed/cxcr4_clean.csv")
  colnames(data)[1:2] <- c("resno", "sub") 
  
  # add resno
  rep <- data.frame()
  for (i in 2:352){
    temp <- as.data.frame(rep(i, 21))
    rep <- rbind(temp, rep)
  }  
  colnames(rep) <- c("resno")
  rep <- rep[order(rep$resno), ]
  rep <- as.data.frame(rep)
  colnames(rep) <- c("resno")
  
  data$resno <- rep$resno
  rm(temp, rep, i)
  
  # add resid
  resid <- data %>% filter(Reads == "WT") %>% select(resno, sub)
  data$resid <- resid$sub[match(unlist(data$resno), resid$resno)]
  rm(resid)
  
  # reorder
  data <- data %>% select(resno, resid, sub, Reads, fitc_1, fitc_2, ab_1, ab_2, alexa_1, alexa_2, cxcl12_1, cxcl12_2)
  data <- data %>% filter(Reads != "WT")
  colnames(data)[4] <- c("reads")
  
  # (1.2) ADD GPCRDB -----------------------------------------------------------
  # gpcrdb
  gpcrdb <- read_csv("01_structure_contacts/data/processed/lookup_pdb_to_gnccn_20200918.csv") %>% 
    select(bw_ngo_ckr, clean_ngo_ckr)
  colnames(gpcrdb) <- c("gn", "resno")
  data <- left_join(data, gpcrdb)
  rm(gpcrdb)
  
  # reorder
  data <- data %>% select(resno, resid, sub, gn, reads, 
                          fitc_1, fitc_2, ab_1, ab_2, alexa_1, alexa_2, cxcl12_1, cxcl12_2)
  
  # write output
  write_csv(data, "09_sat_mut_data/data/processed/cxcr4_clean_gpcrdb.csv")

  
##### 2: MEANS - LINEAR SCALE  #################################################
  
  # (2.1) IMPORT, AVERAGES -----------------------------------------------------
  # import, gather
  data <- read_csv("09_sat_mut_data/data/processed/cxcr4_clean_gpcrdb.csv")
  data <- data %>% gather(rep, value, 6:13) 
  data <- data %>% separate(rep, c("sele", "rep"))
  
  # remove NAs
  data <- data %>% filter(!is.na(value))
  data <- data %>% filter(!is.na(gn))
  
  # transform linear scale
  data <- data %>% mutate(value = (2^(data$value) ))
  
  # mean across RESIDUE
  data <- data %>% group_by(resno, sele) %>%
    mutate(res_mean = mean(value, na.rm = TRUE)) %>% ungroup()
  
  # mean across SPECIFIC SUB
  data <- data %>% group_by(resno, sub, sele) %>%
    mutate(sub_mean = mean(value, na.rm = TRUE)) %>% ungroup()
  
  # re-transform to log2
  data <- data %>%
    mutate(value_log2 = log2(data$value)) %>%
    mutate(sub_mean_log2 = log2(data$sub_mean)) %>%
    mutate(res_mean_log2 = log2(data$res_mean)) # %>%
    # mutate(sub_sd_log2 = log2(data$sub_sd)) %>%
    # mutate(res_sd_log2 = log2(data$res_sd))
    # 
  
  # write output
  write_csv(data, "09_sat_mut_data/output/cxcr4_clean_gpcrdb_means_linear.csv")
  
  
##### 3: MEANS - LOG SCALE  ####################################################
  
  # (2.1) IMPORT, AVERAGES -----------------------------------------------------
  # import, gather
  data <- read_csv("09_sat_mut_data/data/processed/cxcr4_clean_gpcrdb.csv")
  data <- data %>% gather(rep, value, 6:13) 
  data <- data %>% separate(rep, c("sele", "rep"))
  
  # remove NAs
  data <- data %>% filter(!is.na(value))
  data <- data %>% filter(!is.na(gn))
  
  # transform linear scale
  # data <- data %>% mutate(value = (2^(data$value) ))
  
  # mean across RESIDUE
  data <- data %>% group_by(resno, sele) %>%
    mutate(res_mean = mean(value, na.rm = TRUE), res_sd = sd(value, na.rm = TRUE)) %>% ungroup()
  
  # mean across SPECIFIC SUB
  data <- data %>% group_by(resno, sub, sele) %>%
    mutate(sub_mean = mean(value, na.rm = TRUE), sub_sd = sd(value, na.rm = TRUE)) %>% ungroup()
  
  # re-transform to log2
  # data <- data %>%
  #   mutate(value_log2 = log2(data$value)) %>%
  #   mutate(sub_mean_log2 = log2(data$sub_mean)) %>%
  #   mutate(res_mean_log2 = log2(data$res_mean)) # %>%
  # mutate(sub_sd_log2 = log2(data$sub_sd)) %>%
  # mutate(res_sd_log2 = log2(data$res_sd))
  # 
  
  # write output
  write_csv(data, "09_sat_mut_data/output/cxcr4_clean_gpcrdb_means_log.csv")

  