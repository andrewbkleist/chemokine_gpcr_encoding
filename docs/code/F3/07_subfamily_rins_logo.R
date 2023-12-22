source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CHEMOKINE ----------------------------------------------------------------
# define tier2 sequence positions
pos <- c("NTc.Cm3")
    
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
gridExtra::grid.arrange(p1, p2)
rm(pos, seq, group1, group2, p1, p2)

# ggsave(filename = "ck_logo_NTc.Cm3.pdf",
#        plot = last_plot(), path = "output/F3/",
#        width = 3,
#        height = 7)


# (2) CHEMOKINE ----------------------------------------------------------------
# define tier2 sequence positions
pos <- c("gn1x24", "gn6x58")

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
gridExtra::grid.arrange(p1, p2)
rm(pos, seq, group1, group2, p1, p2)

# ggsave(filename = "ckr_logo_1x24_6x58.pdf",
#        plot = last_plot(), path = "output/F3/",
#        width = 3,
#        height = 7)
  
  