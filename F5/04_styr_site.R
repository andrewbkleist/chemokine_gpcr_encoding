# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
require(ggseqlogo)
library(gridExtra)

################################################################################

# (1) sTYR BINDING SITE ------------------------------------------------------
data <- read_csv("01_structure_contacts/output/RIN_residue.csv") %>% 
  filter(file %in% c("2mpm", "2k05", "6fgp")) %>%
  filter(Res2 == "TYS") %>%
  # filter(ResNum2 %in% c(21)) %>%
  filter(Chain1 != Chain2)

# order
order <- c("CX.1","CX.2","CX.3","CX.4","CX.5","cxb1.1","cxb1.2","cxb1.3","cxb1.4","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.9","cxb1.10","cxb1.11","cxb1.12","cxb1.13","cxb1.14","cxb1.15","cxb1.16","cxb1.17","cxb1.18","cxb1.19","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.1","b1b2.2","b1b2.3","b1b2.4","b1b2.5","b1b2.6","b1b2.7","b1b2.8","b1b2.9","b1b2.10","b1b2.11","b1b2.12","b1b2.13","b1b2.14","b1b2.15","b1b2.16","b1b2.17","b1b2.18","b1b2.19","b1b2.20","b1b2.21","b1b2.22","b1b2.23","b1b2.24","b1b2.25","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.1","b2b3.2","b2b3.3","b2b3.4","b2b3.5","b2b3.6","b2b3.7","b2b3.8","b2b3.9","b2b3.10","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","b3h.5","b3h.6","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10")
data$source_gnccn <- factor(data$source_gnccn, levels = order)
data$file <- factor(data$file, levels = c("6fgp","2mpm","2k05"))

data %>%
  ggplot(aes(source_gnccn, rev(as.factor(rev(ResNum2))))) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  facet_grid(file ~ .)

# ggsave(filename = "styr_rin.pdf", 
#        plot = last_plot(), path = "F5/output/",
#        width = 4,
#        height = 6)

# cxb1.6, cxb1.11, B1.1, b3b3.12

gnccn <- data %>% select(source_gnccn) %>% unique()
gnccn <- gnccn$source_gnccn
gnccn <- as.character(gnccn)

# (2) CUSTOM PIE CHARTS --------------------------------------------------------

# import data
data <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")

# coord polar normalizes to the max value, add something with max value 1.0
# note that B3.3 is not a sTyr contact
gnccn <- c("B3.3", "cxb1.6", "cxb1.11", "B1.1", "b2b3.12")

# select and pie chart - CHEMOKINE
data %>% 
  filter(source_gnccn %in% gnccn) %>% 
  select(source_gnccn, all_para_ck) %>% unique() %>%
  ggplot(aes(x = "", y= all_para_ck)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ source_gnccn)

# ggsave(filename = "styr_motif_paralog_cons_pie.pdf", 
#        plot = last_plot(), path = "F5/output/",
#        width = 6,
#        height = 3)

# cxb1.6, cxb1.11, B1.1, b2b3.12


# (3) CUSTOM SEQ LOGO ----------------------------------------------------------

# manually define positions for seq logo
pos <- gnccn

# order
# order <- c("CX.1","CX.2","CX.3","CX.4","CX.5","cxb1.1","cxb1.2","cxb1.3","cxb1.4","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.9","cxb1.10","cxb1.11","cxb1.12","cxb1.13","cxb1.14","cxb1.15","cxb1.16","cxb1.17","cxb1.18","cxb1.19","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.1","b1b2.2","b1b2.3","b1b2.4","b1b2.5","b1b2.6","b1b2.7","b1b2.8","b1b2.9","b1b2.10","b1b2.11","b1b2.12","b1b2.13","b1b2.14","b1b2.15","b1b2.16","b1b2.17","b1b2.18","b1b2.19","b1b2.20","b1b2.21","b1b2.22","b1b2.23","b1b2.24","b1b2.25","B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.1","b2b3.2","b2b3.3","b2b3.4","b2b3.5","b2b3.6","b2b3.7","b2b3.8","b2b3.9","b2b3.10","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","b3h.5","b3h.6","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10")
# data$source_gnccn <- factor(data$source_gnccn, levels = order)

# load alignment df
ck <- read_csv("02_ck_seq/data/processed/ALL_para_df.csv")
ck <- ck[(names(ck) %in% pos)]
print(colnames(ck))
ck <- unite(ck[,1:ncol(ck)], col = seq_string,  sep = "")

p <- ggseqlogo(ck,  method = 'bits' ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
gridExtra::grid.arrange(p)


# ggsave(filename = "styr_motif_paralog_logo.pdf", 
#        plot = last_plot(), path = "F5/output/",
#        width = 3,
#        height = 6)

