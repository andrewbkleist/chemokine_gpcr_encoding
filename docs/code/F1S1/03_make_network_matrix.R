source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import
data <- read_csv("data/network/Supplementary_Table_1.csv") %>%
  dplyr::select(chemokine, gpcr, ligand_type, interaction_strength) %>%
  filter(!(interaction_strength == "ND"))%>%
  unique()

# order positions
order.ck <- as.factor(toupper(c("ccl1","ccl2","ccl3","ccl3l1","ccl4","ccl4l1",
                                "ccl5","ccl7","ccl8","ccl11","ccl13","ccl14",
                                "ccl15","ccl16","ccl17","ccl18","ccl19","ccl20",
                                "ccl21","ccl22","ccl23","ccl24","ccl25","ccl26",
                                "ccl27","ccl28","cxcl1","cxcl2","cxcl3","cxcl4",
                                "cxcl4l1","cxcl5","cxcl6","cxcl7","cxcl8",
                                "cxcl9","cxcl10","cxcl11","cxcl12","cxcl13",
                                "cxcl14","cxcl16","cxcl17","cx3cl1","xcl1",
                                "xcl2")))
data$chemokine <- factor(data$chemokine, levels = rev(order.ck))

order.ckr <- as.factor(toupper(c("ccr1","ccr2","ccr3","ccr4","ccr5","ccr6",
                                 "ccr7","ccr8","ccr9","ccr10","cxcr1","cxcr2",
                                 "cxcr3","cxcr4","cxcr5","cxcr6","ackr1","ackr2",
                                 "ackr3","ackr4","ccrl2","cx3cr1","xcr1")))
data$gpcr <- factor(data$gpcr, levels = (order.ckr))

# change data class
data$interaction_strength <- as.numeric(data$interaction_strength)

data %>%
  filter(!(is.na(gpcr))) %>%
  unique() %>%
  ggplot(aes(gpcr, chemokine, fill = interaction_strength)) +
  geom_tile() +
  scale_fill_gradient(low="grey90", high="black") +
  geom_text(aes(gpcr, chemokine, label = interaction_strength), size = 3) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = "network_matrix.pdf",
       plot = last_plot(), path = "output/F1S/",
       width = 6,
       height = 6)
# note that some pairings have multiple entries - for these, manually selected
# entries with the higher score values and adjusted colors accordingly for
# final figure

# positive pairings
data <- data %>% 
  dplyr::select(chemokine, gpcr, ligand_type, interaction_strength, evidence_grade) %>%
  filter(interaction_strength %in% c(2,3)) %>%
  dplyr::select(chemokine, gpcr) %>% unique()
