# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# import data, select only CC/CXC
rin <- read_csv("05_integrate/output/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo"))

# chemokine
rin %>%
  filter(cc_cxc_lr_ck >= 0.75) %>% # select subfamily-specific
  select(source_gnccn, cc_para_ck, cxc_para_ck, cc_cxc_lr_ck) %>%
  unique() %>%
  ggplot(aes(cc_para_ck, cxc_para_ck)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  # scale_size(range = c(2, 8)) +
  xlim(0,1) +
  ylim(0,1) +
  theme_minimal()

ggsave(filename = "ck_subfamily_cons.pdf", 
       plot = last_plot(), path = "F3S/output/",
       width = 4,
       height = 4)

# receptor
rin %>%
  filter(cc_cxc_lr_ck >= 0.75) %>% # select subfamily-specific
  select(target_gnccn, cc_para_ckr, cxc_para_ckr, cc_cxc_lr_ckr) %>%
  unique() %>%
  ggplot(aes(cc_para_ckr, cxc_para_ckr)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  # scale_size(range = c(2, 8)) +
  xlim(0,1) +
  ylim(0,1) +
  theme_minimal()

ggsave(filename = "ckr_subfamily_cons.pdf", 
       plot = last_plot(), path = "F3S/output/",
       width = 4,
       height = 4)
  