source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) HEATMAP - CHEMOKINE N-TERM -----------------------------------------------
# import, reformat
data <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv")  

# select only relevant info (motif, protein, class, pct ortho)
data <- data %>% select(motif, protein, class, pct_ortho)

# order chemokines
order.ck <- as.factor(tolower(c("ccl1","ccl2", "ccl3", "ccl3l1",
                                "ccl4","ccl4l1","ccl5","ccl7","ccl8",
                                "ccl11","ccl13","ccl14","ccl15",
                                "ccl16","ccl17","ccl18","ccl19",
                                "ccl20","ccl21","ccl22","ccl23",
                                "ccl24","ccl25","ccl26","ccl27",
                                "ccl28","cxcl1","cxcl2","cxcl3",
                                "cxcl4","cxcl4l1","cxcl5","cxcl6",
                                "cxcl7","cxcl8","cxcl9","cxcl10",
                                "cxcl11","cxcl12","cxcl13", "cxcl14", "cxcl16",
                                "cxcl17", "cx3cl1", "xcl1", "xcl2")))
levels(data$protein)
data$protein <- factor(data$protein, levels = order.ck)

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein, -pct_ortho)
order.motif <- unique(data$motif)
levels(data$motif)
data$motif <- factor(data$motif, levels = rev(order.motif))

#(1.2) PLOTTING --------------------------------------------------------------
# MASTER HEATMAP
data %>%
  # filter(motif == "AA") %>%
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumpurple3")+
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ggsave(filename = "motif_heatmap_ck_nt.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 3.5,
#        height = 6)

# (2) HEATMAP - GPCR N-TERM ----------------------------------------------------
# import, reformat
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")  

# select only relevant info (motif, protein, class, pct ortho)
data <- data %>% select(motif, protein, class, pct_ortho)

# order receptors
order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                 "ccr5","ccr6","ccr7","ccr8","ccr9",
                                 "ccr10","cxcr1","cxcr2","cxcr3",
                                 "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                 "ackr1","ackr2","ackr3","ackr4",
                                 "ccrl2")))

levels(data$protein)
data$protein <- factor(data$protein, levels = order.ckr)

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein, -pct_ortho)
order.motif <- unique(data$motif)
levels(data$motif)
data$motif <- factor(data$motif, levels = rev(order.motif))


#(2.2) PLOTTING ----------------------------------------------------------------

# MASTER HEATMAP
# use below to find position of motif; must move "fill" inside aes to work
# data <- data %>% mutate(dyd = case_when(
#   motif == "DYD" ~ "yes",
#   motif != "DYD" ~ "no"
# ))
# data$dyd <- as.factor(data$dyd)

data %>%
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumpurple3") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ggsave(filename = "motif_heatmap_ckr_nt.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 3,
#        height = 6)


# (3) HEATMAP - GPCR ECL2 ------------------------------------------------------
# import, reformat
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv")  

# select only relevant info (motif, protein, class, pct ortho)
data <- data %>% select(motif, protein, class, pct_ortho)

# order receptors
order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                 "ccr5","ccr6","ccr7","ccr8","ccr9",
                                 "ccr10","cxcr1","cxcr2","cxcr3",
                                 "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                 "ackr1","ackr2","ackr3","ackr4",
                                 "ccrl2")))

levels(data$protein)
data$protein <- factor(data$protein, levels = order.ckr)

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein, -pct_ortho)
order.motif <- unique(data$motif)
levels(data$motif)
data$motif <- factor(data$motif, levels = rev(order.motif))


#(3.2) PLOTTING ----------------------------------------------------------------

# MASTER HEATMAP
# use below to find position of motif; must move "fill" inside aes to work
# data <- data %>% mutate(dyd = case_when(
#   motif == "DYD" ~ "yes",
#   motif != "DYD" ~ "no"
# ))
# data$dyd <- as.factor(data$dyd)

data %>%
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumpurple3") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ggsave(filename = "motif_heatmap_ckr_ecl2.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 3,
#        height = 6)

