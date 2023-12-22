source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import, reformat
data <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
  filter(mer == "mer2") %>%
  filter(mask == "none") %>%
  filter(!(protein %in%  c("ccl2","ccl4l1", "ccl5","ccl7","ccl11","ccl13","ccl14",
                         "ccl15","ccl18","ccl24","ccl26", "ccl28",
                         "cxcl11")))

# select only relevant info (motif, protein, class, pct ortho)
data <- data %>% select(motif, protein, class, pct_ortho)

# order 
data$protein <- factor(data$protein, levels = c("ccl2","ccl4l1", "ccl5","ccl7","ccl11","ccl13","ccl14",
                                                "ccl15","ccl18","ccl24","ccl26", "ccl28",
                                                "cxcl11"))

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein, -pct_ortho)
order.motif <- unique(data$motif)
levels(data$motif)
data$motif <- factor(data$motif, levels = rev(order.motif))

data %>%
  filter(pct_ortho >= 0.5) %>%
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumorchid4") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

###########

# import, reformat
data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv") %>%
  filter(mer == "mer3") %>%
  filter(mask == "none") %>%
  filter(protein %in%  c("ccl27", "ccl28"))

# select only relevant info (motif, protein, class, pct ortho)
data <- data %>% select(motif, protein, class, pct_ortho)

# order 
data$protein <- factor(data$protein, levels = c("ccl2","ccl5","ccl7","ccl8","ccl11","ccl13","ccl14","ccl15","ccl18","ccl24","ccl26","ccl27", "ccl28"))

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein, -pct_ortho)
order.motif <- unique(data$motif)
levels(data$motif)
data$motif <- factor(data$motif, levels = rev(order.motif))

data %>%
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumorchid4") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))





# DELETE
data <- read_csv("07_ck_motif/output/CK_MOTIF_FREQUENCY.csv") %>%
  filter(mer == "mer2") %>%
  filter(mask == "none") %>%
  filter(protein %in%  c("ccl8","ccl5",))

