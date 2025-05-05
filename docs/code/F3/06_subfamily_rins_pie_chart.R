source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
# import data
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv")
  
# select and pie chart - CHEMOKINE
rin %>% #  filter(no_pdb >4) %>% 
  filter(source_gnccn %in% c("NTc.Cm1", "NTc.Cm3", "NTc.Cm4", "CX.4",
                             "b1b2.4", "b1b2.6", "b1b2.9", "b1b2.10", "b1b2.16")) %>%
  select(source_gnccn, cc_cxc_lr_ck) %>% unique() %>%
  ggplot(aes(x = "", y= cc_cxc_lr_ck)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ source_gnccn)
  
  
# ggsave(filename = "ck_cc_vs_cxc_pie.pdf",
#        plot = last_plot(), path = "output/F3/",
#        width = 7,
#        height = 3)

# select and pie chart - RECEPTOR
rin %>% #  filter(no_pdb >4) %>% 
  filter(target_gnccn %in% c("1x24", "1x28", "45x51", "ECL2.Cp6", "5x32", "5x36", "6x58", "7x27", "7x34")) %>%
  select(target_gnccn, cc_cxc_lr_ckr) %>% unique() %>%
  ggplot(aes(x = "", y= cc_cxc_lr_ckr)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  ylim(0,1) +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ target_gnccn)

# ggsave(filename = "ckr_cc_vs_cxc_pie.pdf",
#        plot = last_plot(), path = "output/F3/",
#        width = 7,
#        height = 3)
  
  