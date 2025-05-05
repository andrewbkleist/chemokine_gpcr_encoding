source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import data, select only CC/CXC
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0"))

# chemokine
rin %>%
  filter(cc_cxc_lr_ck >= 0.75) %>% # select subfamily-specific
  dplyr::select(source_gnccn, cc_para_ck, cxc_para_ck, cc_cxc_lr_ck) %>%
  unique() %>%
  # filter(source_gnccn %in% c("B2.6", "H.10","cxb1.7")) %>%
  ggplot(aes(cc_para_ck, cxc_para_ck)) +
  geom_point(shape = 21, colour = "white", fill = "firebrick4", size = 8, stroke = 0.5) +
  xlim(0,1) +
  ylim(0,1) +
  theme_minimal()

# ggsave(filename = "ck_cc_by_cxc_subfamily_cons.pdf",
#        plot = last_plot(), path = "output/F3S1/",
#        width = 4,
#        height = 4)

# receptor
rin %>%
  filter(cc_cxc_lr_ckr >= 0.75) %>% # select subfamily-specific
  dplyr::select(target_gnccn, cc_para_ckr, cxc_para_ckr, cc_cxc_lr_ckr) %>%
  unique() %>%
  filter(target_gnccn %in% c("6x58", "1x24")) %>%
  ggplot(aes(cc_para_ckr, cxc_para_ckr)) +
  geom_point(shape = 21, colour = "white", fill = "grey40", size = 8, stroke = 0.5) +
  xlim(0,1) +
  ylim(0,1) +
  theme_minimal()

# ggsave(filename = "ckr_cc_by_cxc_subfamily_cons.pdf",
#        plot = last_plot(), path = "output/F3S1/",
#        width = 4,
#        height = 4)
