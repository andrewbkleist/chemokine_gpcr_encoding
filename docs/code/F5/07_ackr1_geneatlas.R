source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# RECEPTOR
data <- read_csv("data/variant/geneatlas/processed/CKR_GENEATLAS.csv") %>%
  filter(P_value < 1e-8) %>%
  filter(Gene == "ACKR1") %>%
  filter(AA_Consequences == "G42D")

data$Trait_Description <- factor(data$Trait_Description, levels = rev(data$Trait_Description[order(data$P_value)]))

data %>%
  ggplot(aes(-log10(P_value), Trait_Description)) +
  geom_bar(stat = "identity") +
  theme_minimal()

# ggsave(filename = "ackr1_traits.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 7,
#        height = 3)
