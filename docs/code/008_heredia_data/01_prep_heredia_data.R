# Imports raw data from Heredia, et al (PMID 29678950) and prepares data for
# subsequent analysis. Downloaded data files manually from paper URL
# (https://journals.aai.org/jimmunol/article/200/11/3825/106401/Mapping-Interaction-Sites-on-Human-Chemokine)
# "saved "cleaned" data by removing annotations, leaving only raw data, saved 
# as cxcr4_clean.csv. 
#
# Emailed Erik Procko (dataset corr. author) to clarify "CXCL12 binding 
# conservation score" and "Binding conservation score" meaning from Fig 3C 
# (response email dated 9/8/2021):
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
# Procko suggested averaging over all mutations:
# "Whenever we average replicate data sets we treat it as combining the replicates 
# to give a single pooled sample, so averaging is with the raw ratios.  We only 
#  average the log2 ratios when calculating a conservation score, which is an 
# abstract term to communicate whether a position is mutationally tolerant or not.
# ...I very much view the CCR5/CXCR4 data as more predictive and qualitative"
#
# Note regarding "NA" values from paper:
# "NA" values reflect missing values: "Missing mutations (<10 reads in the naive
# library) are black." (Fig 1)

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) IMPORT, ADD RESNO, RESID -------------------------------------------------
# import data
data <- read_csv("data/mutagenesis_heredia/raw/cxcr4_clean.csv")
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
resid <- data %>% filter(Reads == "WT") %>% dplyr::select(resno, sub)
data$resid <- resid$sub[match(unlist(data$resno), resid$resno)]
rm(resid)

# reorder
data <- data %>% dplyr::select(resno, resid, sub, Reads, fitc_1, fitc_2, ab_1, ab_2, alexa_1, alexa_2, cxcl12_1, cxcl12_2)
data <- data %>% filter(Reads != "WT")
colnames(data)[4] <- c("reads")

# (2) ADD GPCRDB ---------------------------------------------------------------
# gpcrdb
gpcrdb <- read_csv("data/lookup/lookup_pdb_to_gnccn_20230924.csv") %>% 
  dplyr::select(bw_ngo_ckr, clean_ngo_ckr)
colnames(gpcrdb) <- c("gn", "resno")
data <- left_join(data, gpcrdb)
rm(gpcrdb)

# reorder
data <- data %>% dplyr::select(resno, resid, sub, gn, reads, 
                        fitc_1, fitc_2, ab_1, ab_2, alexa_1, alexa_2, cxcl12_1, cxcl12_2)

# write output
# write_csv(data, "data/mutagenesis_heredia/raw/cxcr4_clean_gpcrdb.csv") # WRITTEN 20231106

# (3) TAKE MEANS - LINEAR SCALE ------------------------------------------------
# import, gather
data <- read_csv("data/mutagenesis_heredia/raw/cxcr4_clean_gpcrdb.csv")
data <- data %>% gather(rep, value, 6:13) 
data <- data %>% separate(rep, c("sele", "rep"))

# remove NAs
data <- data %>% filter(!is.na(value))
data <- data %>% filter(!is.na(gn))

# transform linear scale
data <- data %>% mutate(value = (2^(data$value) ))

# mean across RESIDUE
data <- data %>% group_by(resno, sele) %>%
  dplyr::mutate(res_mean = mean(value, na.rm = TRUE)) %>% ungroup()

# mean across SPECIFIC SUB
data <- data %>% group_by(resno, sub, sele) %>%
  dplyr::mutate(sub_mean = mean(value, na.rm = TRUE)) %>% ungroup()

# re-transform to log2
data <- data %>%
  dplyr::mutate(value_log2 = log2(data$value)) %>%
  dplyr::mutate(sub_mean_log2 = log2(data$sub_mean)) %>%
  dplyr::mutate(res_mean_log2 = log2(data$res_mean)) 

# write output
# write_csv(data, "data/mutagenesis_heredia/processed/cxcr4_clean_gpcrdb_means_linear.csv") # WRITTEN 20231106

# (4) TAKE MEANS - LOG SCALE ---------------------------------------------------
# import, gather
data <- read_csv("data/mutagenesis_heredia/raw/cxcr4_clean_gpcrdb.csv")
data <- data %>% gather(rep, value, 6:13) 
data <- data %>% separate(rep, c("sele", "rep"))

# remove NAs
data <- data %>% filter(!is.na(value))
data <- data %>% filter(!is.na(gn))

# mean across RESIDUE
data <- data %>% group_by(resno, sele) %>%
  dplyr::mutate(res_mean = mean(value, na.rm = TRUE), res_sd = sd(value, na.rm = TRUE)) %>% ungroup()

# mean across SPECIFIC SUB
data <- data %>% group_by(resno, sub, sele) %>%
  dplyr::mutate(sub_mean = mean(value, na.rm = TRUE), sub_sd = sd(value, na.rm = TRUE)) %>% ungroup()

# write output
# write_csv(data, "data/mutagenesis_heredia/processed/cxcr4_clean_gpcrdb_means_log.csv") # WRITTEN 20231106
