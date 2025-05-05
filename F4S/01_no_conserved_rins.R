# Name:     ***
# Updated:  ***
# Author:   Andrew Kleist

# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv")
rin <- rin %>% 
  filter(file != "6meo") %>%
  # filter(source_gnccn == "cxb1.1" & target_gnccn == "7x27") %>%
  count(source_gnccn, target_gnccn) %>% count(n)

rin %>%
  ggplot(aes(rev(n), nn)) +
  geom_bar(stat = "identity") +
  theme_minimal()

ggsave(filename = "no_complexes_bar.pdf", 
       plot = last_plot(), path = "F4S/output/",
       width = 4,
       height = 4)
