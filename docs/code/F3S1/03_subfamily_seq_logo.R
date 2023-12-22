source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CHEMOKINE ----------------------------------------------------------------
pos <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0")) %>%
  filter(cc_cxc_lr_ck >= 0.75) %>% # select subfamily-specific
  dplyr::select(source_gnccn) %>% unique()
pos <- pos$source_gnccn

# load alignment df
seq <- read_csv("data/sequence/chemokine/alignment_csv/ALL_para_df.csv")

# make CC/CXC and unique columns
group1 <- seq %>% filter(class == "cc")
group1 <- group1[(names(group1) %in% pos)]
group1 <- unite(group1[,1:ncol(group1)], col = seq_string,  sep = "")

group2 <- seq %>% filter(class == "cxc")
group2 <- group2[(names(group2) %in% pos)]
group2 <- unite(group2[,1:ncol(group2)], col = seq_string,  sep = "")

# plot
p1 <- ggseqlogo(group1,  method = 'bits' ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- ggseqlogo(group2,  method = 'bits' ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- gridExtra::grid.arrange(p1, p2)
rm(pos, seq, group1, group2, p1, p2)

# ggsave(filename = "ck_subfamily_logo.pdf",
#        plot = p, path = "output/F3S1/",
#        width = 8,
#        height = 4)
rm(p)

# (2) RECEPTOR -----------------------------------------------------------------
pos <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file %in% c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", 
                     "7xa3", "7f1t", "6wwz", "6lfo", "ngo", "8ic0")) %>%
  filter(cc_cxc_lr_ckr >= 0.75) %>% # select subfamily-specific
  dplyr::select(target_gnccn) %>% unique()
pos$gn <- c("gn")
pos <- pos %>% unite(col = target_gnccn, c(gn, target_gnccn), sep = "")
pos <- pos$target_gnccn

# load alignment df
seq <- read_csv("data/sequence/gpcr/alignment_csv/ALL_classa_df.csv")

# make CC/CXC and unique columns
group1 <- seq %>% filter(class == "cc")
group1 <- group1[(names(group1) %in% pos)]
group1 <- unite(group1[,1:ncol(group1)], col = seq_string,  sep = "")

group2 <- seq %>% filter(class == "cxc")
group2 <- group2[(names(group2) %in% pos)]
group2 <- unite(group2[,1:ncol(group2)], col = seq_string,  sep = "")

# plot
p1 <- ggseqlogo(group1,  method = 'bits' ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p2 <- ggseqlogo(group2,  method = 'bits' ) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
p <- gridExtra::grid.arrange(p1, p2)
p
rm(pos, seq, group1, group2, p1, p2)

# ggsave(filename = "ckr_subfamily_logo.pdf",
#        plot = p, path = "output/F3S1/",
#        width = 8,
#        height = 4)
rm(p)