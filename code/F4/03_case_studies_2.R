source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) IDENTIFY CONSENSUS CXCL8 CONTACTS ----------------------------------------
# 1a: import rin comparison
pdbs <- c("8ic0", "6lfo") # CXCL8:CXCR2 / CXCR8:CXCR1
rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  # filter(type %in% c("a_not_b","b_not_a")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  mutate(str_unstr_ck = case_when(
    grepl("NTc",source_gnccn) ~ "unstr",
    !grepl("NTc",source_gnccn) ~ "str")) %>%
  mutate(str_unstr_ckr = case_when(
    grepl("NTr",target_gnccn) ~ "unstr",
    !grepl("NTr",target_gnccn) ~ "str"))
rm(pdbs)

# 1b: generate seq comparison
seq <- read_csv("data/sequence/gpcr/alignment_csv/ALL_classa_df.csv") %>%
  filter(protein %in% c("cxcr1","cxcr2")) %>%
  dplyr::select(-class, -seq)
seq <- seq %>%
  pivot_longer(cols = 2:ncol(seq), names_to = "target_gnccn", values_to = "resid")
seq <- seq %>% pivot_wider(names_from = protein, values_from = resid)  
seq <- seq %>%
  mutate(identity = case_when(
    cxcr1 == cxcr2 ~ 1,
    cxcr1!= cxcr2 ~ 0
  ))
seq <- seq %>% separate(col = target_gnccn, into = c("a", "target_gnccn"), sep = "gn")
seq <- seq %>% dplyr::select(-a)
rin <- left_join(rin, seq)
rm(seq)

# 1c: now add ortholog conservation info
lookup <- read_csv("data/integrated/CK_CONS_CCCXC_SNP_CAN.csv") %>%
  filter(protein == "CXCL8") %>% 
  dplyr::select(ccn, ortho_cons)
colnames(lookup) <- c("source_gnccn", "ortho_cons_cxcl8")  
rin$ortho_cons_cxcl8 <- lookup$ortho_cons_cxcl8[match(unlist(rin$source_gnccn), lookup$source_gnccn)]
rm(lookup)

lookup <- read_csv("data/integrated/CKR_CONS_CCCXC_SNP_CAN.csv") %>%
  filter(protein == "CXCR1") %>% 
  dplyr::select(gn, ortho_cons) %>%
  mutate(gn = case_when(
    gn == "gnECL3.5" ~ "gnECL3.Cm2",
    gn != "gnECL3.5" ~ gn
    
  ))
colnames(lookup) <- c("target_gnccn", "ortho_cons_cxcr1")  
lookup <- lookup %>% separate(col = target_gnccn, into = c("gn", "target_gnccn"), sep  = "gn") %>% dplyr::select(-gn)
rin$ortho_cons_cxcr1 <- lookup$ortho_cons_cxcr1[match(unlist(rin$target_gnccn), lookup$target_gnccn)]
rm(lookup)

lookup <- read_csv("data/integrated/CKR_CONS_CCCXC_SNP_CAN.csv") %>%
  filter(protein == "CXCR2") %>% 
  dplyr::select(gn, ortho_cons)  %>%
  mutate(gn = case_when(
    gn == "gnECL3.5" ~ "gnECL3.Cm2",
    gn != "gnECL3.5" ~ gn
    ))
colnames(lookup) <- c("target_gnccn", "ortho_cons_cxcr2")  
lookup <- lookup %>% separate(col = target_gnccn, into = c("gn", "target_gnccn"), sep  = "gn") %>% dplyr::select(-gn)
rin$ortho_cons_cxcr2 <- lookup$ortho_cons_cxcr2[match(unlist(rin$target_gnccn), lookup$target_gnccn)]
rm(lookup)

# 1d: identify shared contacts that use identical residues
# note of the two for which there are NAs, ECL2.Cm9 are different 
# (Y178 for CXCR1 and V187 for CXCR2) but ECL3.Cm2 are shared (E275 for CXCR1
# and E284 for CXCR2)
rin <- rin %>%
  mutate(identity = case_when(
    target_gnccn == "ECL3.Cm2" ~ 1,
    target_gnccn != "ECL3.Cm2" ~ identity
  )) %>%
  mutate(identity = case_when(
    target_gnccn == "ECL2.Cm9" ~ 0,
    target_gnccn != "ECL2.Cm9" ~ identity 
  ))
rin <- rin %>% filter(identity == 1) %>% 
  filter(type == "a_and_b") %>%
  filter(ortho_cons_cxcl8 >= 0.5) %>%
  filter(ortho_cons_cxcr1 >= 0.5) %>%
  filter(ortho_cons_cxcr2 >= 0.5)
paste0("CXCL8 determinants are as follows:")
rin

#--
temp <- rin %>% dplyr::count(source_gnccn)
temp <- rin %>% dplyr::count(target_gnccn)
#--

# write conect
conect <- rin %>% dplyr::select(source_gnccn, target_gnccn)
WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  rin,
                  "6lfo",
                  "data/structure/raw/pdb_files/6lfo_clean.pdb",
                  "output/F4/xx_conect/cxcl8_determinants.csv")

# (2) IDENTIFY CONSENSUS CCR5 CONTACTS -----------------------------------------
# 2a: import rin comparison
pdbs <- c("7f1t", "zheng")
rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  mutate(str_unstr_ck = case_when(
    grepl("NTc",source_gnccn) ~ "unstr",
    !grepl("NTc",source_gnccn) ~ "str")) %>%
  mutate(str_unstr_ckr = case_when(
    grepl("NTr",target_gnccn) ~ "unstr",
    !grepl("NTr",target_gnccn) ~ "str"))
rm(pdbs)

# 2b: generate seq comparison
seq <- read_csv("data/sequence/chemokine/alignment_csv/ALL_para_df.csv") %>%
  filter(protein %in% c("ccl3","ccl5")) %>%
  dplyr::select(-class, -seq)
seq <- seq %>%
  pivot_longer(cols = 2:ncol(seq), names_to = "source_gnccn", values_to = "resid")
seq <- seq %>% pivot_wider(names_from = protein, values_from = resid)  
seq <- seq %>%
  mutate(identity = case_when(
    ccl3 == ccl5 ~ 1,
    ccl3!= ccl5 ~ 0
  ))
rin <- left_join(rin, seq)
rm(seq)
# TOTAL:          127
# CCL3-CCR5 only: 40 (31%)
# CCL5-CCR5 only: 41 (32%)
# shared:         46 (36%)
# CCR5:           31 (25%)

# 2c: 
lookup <- read_csv("data/integrated/CK_CONS_CCCXC_SNP_CAN.csv") %>%
  filter(protein == "CCL3") %>% 
  dplyr::select(ccn, ortho_cons)
colnames(lookup) <- c("source_gnccn", "ortho_cons_ccl3")  
rin$ortho_cons_ccl3 <- lookup$ortho_cons_ccl3[match(unlist(rin$source_gnccn), lookup$source_gnccn)]
rm(lookup)

lookup <- read_csv("data/integrated/CK_CONS_CCCXC_SNP_CAN.csv") %>%
  filter(protein == "CCL5") %>% 
  dplyr::select(ccn, ortho_cons)
colnames(lookup) <- c("source_gnccn", "ortho_cons_ccl5")  
rin$ortho_cons_ccl5 <- lookup$ortho_cons_ccl5[match(unlist(rin$source_gnccn), lookup$source_gnccn)]
rm(lookup)

lookup <- read_csv("data/integrated/CKR_CONS_CCCXC_SNP_CAN.csv") %>%
  filter(protein == "CCR5") %>% 
  dplyr::select(gn, ortho_cons) %>%
  mutate(gn = case_when(
    gn == "gnECL2.6" ~ "gnECL2.Cm8",
    gn != "gnECL2.6" ~ gn
  )) %>%
  mutate(gn = case_when(
    gn == "gnECL2.7" ~ "gnECL2.Cm7",
    gn != "gnECL2.7" ~ gn
  )) %>%
  mutate(gn = case_when(
    gn == "gnECL2.8" ~ "gnECL2.Cm6",
    gn != "gnECL2.8" ~ gn
  )) %>%
  mutate(gn = case_when(
    gn == "gnECL3.5" ~ "gnECL3.Cm2",
    gn != "gnECL3.5" ~ gn
  )) 
colnames(lookup) <- c("target_gnccn", "ortho_cons_ccr5")  
lookup <- lookup %>% separate(col = target_gnccn, into = c("gn", "target_gnccn"), sep  = "gn") %>% dplyr::select(-gn)
rin$ortho_cons_ccr5 <- lookup$ortho_cons_ccr5[match(unlist(rin$target_gnccn), lookup$target_gnccn)]
rm(lookup)

# 2d: identify shared contacts using same residues with ortho cons
rin <- rin %>% filter(identity == 1) %>% 
  filter(type == "a_and_b") %>%
  filter(ortho_cons_ccl3 >= 0.5) %>%
  filter(ortho_cons_ccl5 >= 0.5) %>%
  filter(ortho_cons_ccr5 >= 0.5) # %>%
  # filter(!(source_gnccn %in% c("CX.1","b1b2.12","B3.3"))) %>%
  # filter(!(target_gnccn %in% c("1x22", "NTr.Cm1")))
paste0("CCR5 determinants are as follows:")
rin

#--
temp <- rin %>% dplyr::count(source_gnccn)
temp <- rin %>% dplyr::count(target_gnccn)
#--


# write conect
conect <- rin %>% dplyr::select(source_gnccn, target_gnccn)
WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  conect,
                  "7f1t",
                  "data/structure/raw/pdb_files/7f1t_clean.pdb",
                  "output/F4/xx_conect/ccr5_determinants.csv")





#--

##### (1) IDENTIFY  SSE PAIRINGS -----------------------------------------------
# pdbs <- c("7f1t", "zheng") # CCL3:CCR5 / CCL5:CCR5
# pdbs <- c("7f1t", "6wwz") # CCL3:CCR5 / CCL20:CCR6
# pdbs <- c("5uiw", "7o7f") # CCL5 [5P7]:CCR5 / CCL5 [6P4]:CCR5
# pdbs <- c("8ic0", "6lfo") # CXCL8:CXCR2 / CXCR8:CXCR1
pdbs <- c("8ic0", "ngo") # CXCL8:CXCR2 / CXCR8:CXCR1

# pdbs <- c("7f1t", "7vl9") # CCL3:CCR5 / CCL15:CCR1
# pdbs <- c("5uiw", "7o7f", "7f1r", "zheng",
#           "7vl9", "7xa3", "7f1t", "6wwz",
#           "8ic0", "6lfo", "ngo",
#           "7sk3", "7xbx")
# pdbs <- c("5uiw", "7o7f", "7f1r", "zheng",
#           "7vl9", "7xa3", "7f1t") # intra-subfamily, in-network
# pdbs <- c("8ic0", "6lfo", "ngo") # intra-subfamily, in-network


a <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  # filter(type %in% c("a_not_b","b_not_a")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%

  mutate(str_unstr_ck = case_when(
    grepl("NTc",source_gnccn) ~ "unstr",
    !grepl("NTc",source_gnccn) ~ "str")) %>%
  mutate(str_unstr_ckr = case_when(
    grepl("NTr",target_gnccn) ~ "unstr",
    !grepl("NTr",target_gnccn) ~ "str"))

a <- a %>% dplyr::count(type, dom1, dom2)

# order domains
order.ck.dom <- as.factor(unique(c("NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","CX","CX","CX","CX","CX","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","B1","B1","B1","B1","B1","B1","B1","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","B2","B2","B2","B2","B2","B2","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","B3","B3","B3","B3","b3h1","b3h1","b3h1","b3h1","b3h1","b3h1","H","H","H","H","H","H","H","H","H","H", "CT")))
a$dom1 <- factor(a$dom1, levels = rev(order.ck.dom))
order.ckr.dom <- as.factor(unique(c("NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","ICL1","ICL1","ICL1","ICL1","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","ECL1","ECL1","ECL1","ECL1","ECL1","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","ICL3","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","ECL3","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","CT")))
a$dom2 <- factor(a$dom2, levels = order.ckr.dom) 

# PLOT SSE MATRIX
a %>%
  ggplot(aes(dom2, dom1, fill = n)) +
  geom_tile() +
  scale_fill_gradient(low="grey90", high="black") +
  geom_text(aes(dom2, dom1, label = n), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(. ~ type)


##### (2) IDENTIFY SPECIFIC SSE ON CK/CKR -------------------------------------- 
# pdbs <- c("7f1t", "zheng") # CCL3:CCR5 / CCL5:CCR5
data <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  # filter(type %in% c("a_not_b","b_not_a")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  
  mutate(str_unstr_ck = case_when(
    grepl("NTc",source_gnccn) ~ "unstr",
    !grepl("NTc",source_gnccn) ~ "str")) %>%
  mutate(str_unstr_ckr = case_when(
    grepl("NTr",target_gnccn) ~ "unstr",
    !grepl("NTr",target_gnccn) ~ "str"))

a <- data %>% dplyr::count(type, dom1)
b <- data %>% dplyr::count(type, dom2)

a %>%
  ggplot(aes(dom1, n)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  facet_grid(. ~ type)

b %>%
  ggplot(aes(dom2, n)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  facet_grid(. ~ type)

##### (3) WRITE CONECT ---------------------------------------------------------
pdbs <- c("8ic0", "6lfo") # CXCL8:CXCR2 / CXCR8:CXCR1
rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("a_not_b")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  dplyr::select(source_gnccn, target_gnccn)

WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                            rin,
                            "6lfo",
                            "data/structure/raw/pdb_files/6lfo_clean.pdb",
                            "output/F4/xx_conect/6lfo_not_8ic0.csv")

rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("b_not_a")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  dplyr::select(source_gnccn, target_gnccn)

WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  rin,
                  "8ic0",
                  "data/structure/raw/pdb_files/8ic0_clean.pdb",
                  "output/F4/xx_conect/8ic0_not_6lfo.csv")

rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("a_and_b")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  dplyr::select(source_gnccn, target_gnccn)

WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  rin,
                  "8ic0",
                  "data/structure/raw/pdb_files/8ic0_clean.pdb",
                  "output/F4/xx_conect/8ic0_and_6lfo.csv")

##### (4) COMPAREE SEQS --------------------------------------------------------
data <- read_csv("data/sequence/gpcr/alignment_csv/ALL_classa_df.csv") %>%
  filter(protein %in% c("cxcr1","cxcr2")) %>%
  dplyr::select(-class, -seq)

data <- data %>%
  pivot_longer(cols = 2:ncol(data), names_to = "position", values_to = "resid")
data <- data %>% pivot_wider(names_from = protein, values_from = resid)  
data <- data %>%
  mutate(idenity = case_when(
    cxcr1 == cxcr2 ~ 1,
    cxcr1!= cxcr2 ~ 0
  ))


##### (5) WRITE CONECT CCR5 ---------------------------------------------------------
pdbs <- c("7f1t", "zheng")
rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("a_not_b")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  dplyr::select(source_gnccn, target_gnccn)

WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  rin,
                  "7f1t",
                  "data/structure/raw/pdb_files/7f1t_clean.pdb",
                  "output/F4/xx_conect/7f1t_not_zheng.csv")

rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("b_not_a")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  dplyr::select(source_gnccn, target_gnccn)

WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  rin,
                  "zheng",
                  "data/structure/raw/pdb_files/zheng_model_clean.pdb",
                  "output/F4/xx_conect/zheng_not_7f1t.csv")

rin <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("a_and_b")) %>%
  filter(pdb1 %in% c(pdbs)) %>% 
  filter(pdb2 %in% c(pdbs)) %>%
  dplyr::select(source_gnccn, target_gnccn)

WriteCONECTcustom("data/structure/processed/RIN_residue.csv",
                  rin,
                  "7f1t",
                  "data/structure/raw/pdb_files/7f1t_clean.pdb",
                  "output/F4/xx_conect/7f1t_and_zheng.csv")
