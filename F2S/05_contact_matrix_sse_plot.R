# Description:
# Calculates mean number of residue-reside contacts across all SSE pairs among
# chemokine-GPCR complex structures. First calculates number of pairwise domain 
# interactions individually, (ie per complex) THEN get percentages (per complex) 
# to normalize percent of that domain-domain interaction comprising total 
# interface, then get the averages of the normalized pairwise contatcs. 
# NOTE: Because the Ngo model encompasses the entire interaction, all Nterm 
# residues of the receptor from NTr.Cm13 and prior were excluded so that the 
# receptor N-terminus approximates the length of the other complexes. NTr.Cm12 
# is the longest length of another receptor, hence exclusion of NTr.Cm13

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# import contacts
rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>%
  filter(!( target_gnccn %in% c("NTr.Cm27", "NTr.Cm26", "NTr.Cm25", "NTr.Cm24",
                                "NTr.Cm23","NTr.Cm22","NTr.Cm21","NTr.Cm20",
                                "NTr.Cm19","NTr.Cm18","NTr.Cm17","NTr.Cm16",
                                "NTr.Cm15","NTr.Cm14","NTr.Cm13"))) %>%
  filter(file != "6meo") # remove CCR5-gp120 complex

# count number SSE-to-SSE interactions per file
dom <- rin %>% select(dom1, dom2, file) %>% count(dom1, dom2, file)

# get ALL paired SSE-to-SSE counts for ALL complexes
# (ie add zeros to SSEs which do not have representative contact)
dom <- dom %>% unite(col = dom, dom1:dom2, sep = "_")
dom <- dom %>% pivot_wider(names_from = file, values_from = "n")
dom[is.na(dom)] <- 0
dom <- dom %>% pivot_longer(cols = 2:ncol(dom), names_to = "file", values_to = "n" )
dom <- dom %>% separate(col = dom, into = c("dom1", "dom2"), sep = "_")

# count totals for each group
dom <- dom %>% group_by(file) %>% mutate(sum = sum(n)) %>% ungroup()

# get mean no. contacts per SSE pair across all structures
dom <- dom %>% group_by(dom1, dom2) %>% summarize(mean = mean(n), sd = sd(n)) %>% ungroup()

# order domains
order.ck.dom <- as.factor(unique(c("NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","NTc","CX","CX","CX","CX","CX","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","cxb1","B1","B1","B1","B1","B1","B1","B1","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","b1b2","B2","B2","B2","B2","B2","B2","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","b2b3","B3","B3","B3","B3","b3h1","b3h1","b3h1","b3h1","b3h1","b3h1","H","H","H","H","H","H","H","H","H","H", "CT")))
dom$dom1 <- factor(dom$dom1, levels = rev(order.ck.dom))
order.ckr.dom <- as.factor(unique(c("NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","NTr","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","TM1","ICL1","ICL1","ICL1","ICL1","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","TM2","ECL1","ECL1","ECL1","ECL1","ECL1","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","TM3","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","ICL2","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","TM4","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","ECL2","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","TM5","ICL3","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","TM6","ECL3","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","TM7","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","H8","CT")))
dom$dom2 <- factor(dom$dom2, levels = order.ckr.dom)

# round means
dom <- dom %>% mutate(mean = round(mean, digits = 1))

# PLOT SSE MATRIX
dom %>%
  ggplot(aes(dom2, dom1, fill = mean)) +
  geom_tile() +
  scale_fill_gradient(low="grey90", high="black") +
  geom_text(aes(dom2, dom1, label = mean), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "contact_matrix_sse_plot.pdf", 
 plot = last_plot(), path = "F2/output/",
 width = 5,
 height = 4)
