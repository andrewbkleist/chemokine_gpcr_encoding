source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# select and pie chart - CHEMOKINE
read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(source_gnccn %in% c("B2.3", "B3.3", "CX.5", "CX.1", "b1b2.12")) %>%
  dplyr::select(source_gnccn, all_para_ck) %>% 
  unique() %>%
  ggplot(aes(x = "", y= all_para_ck)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  ylim(0,1) + # **THIS IS ESSENTIAL TO INCLUDE FOR ALL PIE CHARTS**
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ source_gnccn)

# ggsave(filename = "conserved_rins_pie_ck.pdf",
#        plot = last_plot(), path = "output/F2/",
#        width = 6,
#        height = 3)

# For conservation scores for plot:
read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(source_gnccn %in% c("B2.3", "B3.3", "CX.5", "CX.1", "b1b2.12")) %>%
  dplyr::select(source_gnccn, all_para_ck) %>% 
  unique()

# select and pie chart - RECEPTOR
read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(target_gnccn %in% c("NTr.Cm1", "1x22", "7x24")) %>%
  dplyr::select(target_gnccn, all_non_ackr_para_ckr) %>% 
  unique() %>%
  ggplot(aes(x = "", y= all_non_ackr_para_ckr)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  ylim(0,1) + # **THIS IS ESSENTIAL TO INCLUDE FOR ALL PIE CHARTS**
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ target_gnccn)

# ggsave(filename = "conserved_rins_pie_ckr.pdf",
#        plot = last_plot(), path = "output/F2/",
#        width = 3,
#        height = 3)

# For conservation scores for plot:
read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(target_gnccn %in% c("NTr.Cm1", "1x22", "7x24")) %>%
  dplyr::select(target_gnccn, all_non_ackr_para_ckr) %>% 
  unique()
