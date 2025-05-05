# Name:     ***
# Updated:  ***
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

ck <- read_csv("02_ck_seq/output/CK_CONSERVATION.csv") %>%
  filter(ccn %in% c("cxb1.1", "cxb1.2", "B3.1"))
ck$ccn <- factor(ck$ccn, levels = c("cxb1.1", "cxb1.2", "B3.1"))

ck %>%
  select(protein, ccn, ortho_cons) %>% 
  unique() %>%
  ggplot(aes(ccn, ortho_cons)) +
  # geom_boxplot() +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2) +
  theme_minimal()

ggsave(filename = "ck_ortholog_cons.pdf", 
       plot = last_plot(), path = "F4S/output/",
       width = 3,
       height = 4)

ckr <- read_csv("03_ckr_seq/output/CKR_CONSERVATION.csv") %>%
  filter(gn %in% c("gnNTr.Cm2", "gnNTr.Cm3", "gn7x27"))
ckr$gn <- factor(ckr$gn, levels = c("gnNTr.Cm2", "gnNTr.Cm3", "gn7x27"))

ckr %>%
  select(protein, gn, ortho_cons) %>% 
  unique() %>%
  ggplot(aes(gn, ortho_cons)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2) +
  theme_minimal()

ggsave(filename = "ckr_ortholog_cons.pdf", 
       plot = last_plot(), path = "F4S/output/",
       width = 3,
       height = 4)

