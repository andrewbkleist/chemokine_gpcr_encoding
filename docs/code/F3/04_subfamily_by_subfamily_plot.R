source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
# import data, select only CC/CXC
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", "7xa3", 
                       "7f1t", "6wwz", "6lfo", "ngo", "8ic0"))
  
# plot (with sizing by no_pdb)
rin %>%
  ggplot(aes(cc_cxc_lr_ckr, cc_cxc_lr_ck))  +
  geom_point(shape = 21, colour = "black", fill = "white", size = 5, stroke = 0.5) +
  xlim(0.5,1) +
  ylim(0.5,1) +
  theme_minimal()

# note that missing values come from contacts made by positions in 5uiw and 
# 7o7f N-termini, which are non-native (i.e. are modified relative to CCL5) and
# ECL2.Cm* positions, which are aligned in such a way that they are non-
# structurally equivalent

# ggsave(filename = "subfamily_by_subfamily.pdf",
#        plot = last_plot(), path = "output/F3/",
#        width = 4.3,
#        height = 4)
  