# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

#################################################################################

data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv") %>%
  filter(protein %in% c("ccl2", "ccl3", "ccl4", "ccl5","ccl7","ccl8","ccl11","ccl13","ccl14","ccl15","ccl18","ccl24","ccl26", "ccl28", "cxcl12"))

# CCR3: c("ccl2","ccl5","ccl7","ccl8","ccl11","ccl13","ccl14","ccl15","ccl18","ccl24","ccl26", "ccl28")
# CCR5: 3,4,5,7,8,11,13,14,26 (NOT 15, 18, 24, 28),

# select only relevant info (motif, protein, class, pct ortho)
data <- data %>% select(motif, protein, class, pct_ortho)

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein, -pct_ortho)
order.motif <- unique(data$motif)
levels(data$motif)
data$motif <- factor(data$motif, levels = rev(order.motif))

# order chemokines
data$protein <- factor(data$protein, levels = (c("ccl2", "ccl3", "ccl4", "ccl5","ccl7","ccl8","ccl11","ccl13","ccl14","ccl15","ccl18","ccl24","ccl26", "ccl28", "cxcl12")))


data %>%
  filter(motif %in% c("PxT", "IP")) %>%
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumorchid4") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))


ggsave(filename = "vmipii_motifs.pdf", 
       plot = last_plot(), path = "F6S/output/",
       width = 4,
       height = 3)


# (1.0) HOW MANY SLiMs? --------------------------------------------------------
# count how many unique fragments there are
data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv")