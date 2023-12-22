source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import data
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv")

# select and pie chart - CHEMOKINE, ORTHOLOG
rin %>% 
  filter(source_gnccn %in% c("cxb1.1")) %>%
  filter(file %in% c("7f1t", "5uiw", "7xa3", "6lfo", "ngo", "7xbx")) %>%
  select(file, ortho_cons_ck) %>% unique() %>%
  ggplot(aes(x = "", y= ortho_cons_ck)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ file)


# ggsave(filename = "cxb1.1_orthol_pie.pdf",
#        plot = last_plot(), path = "output/F4/",
#        width = 6,
#        height = 3)

# select and pie chart - CHEMOKINE, PARALOG
rin %>% 
  filter(source_gnccn %in% c("cxb1.1")) %>%
  select(source_gnccn, all_para_ck) %>% unique() %>%
  ggplot(aes(x = "", y= all_para_ck)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ source_gnccn)

# ggsave(filename = "cxb1.1_paralog_pie.pdf",
#        plot = last_plot(), path = "output/F4/",
#        width = 3,
#        height = 3)

# select and pie chart - RECEPTOR, ORTHOLOG
rin %>% 
  filter(target_gnccn %in% c("7x27")) %>%
  filter(file %in% c("7f1t", "5uiw", "7xa3", "6lfo", "ngo", "7xbx")) %>%
  select(file, ortho_cons_ckr) %>% unique() %>%
  ggplot(aes(x = "", y= ortho_cons_ckr)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ file)

# ggsave(filename = "7x27_ortholog_pie.pdf",
#        plot = last_plot(), path = "output/F4/",
#        width = 6,
#        height = 3)

# select and pie chart - RECEPTOR, PARALOG
rin %>% 
  filter(target_gnccn %in% c("7x27")) %>%
  # filter(file %in% c("7f1t", "5uiw", "7xa3", "6lfo", "ngo", "7xbx")) %>%
  select(target_gnccn, all_non_ackr_para_ckr) %>% unique() %>%
  ggplot(aes(x = "", y= all_non_ackr_para_ckr)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  ylim(0,1) +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) +
  facet_grid(. ~ target_gnccn)

# ggsave(filename = "7x27_paralog_pie.pdf",
#        plot = last_plot(), path = "output/F4/",
#        width = 3,
#        height = 3)
