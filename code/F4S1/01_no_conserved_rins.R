source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo") %>%
  dplyr::count(source_gnccn, target_gnccn) %>%
  dplyr::count(n)

rin %>%
  ggplot(aes(n, nn)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(breaks = seq(0, 12, by = 2)) +
  theme_minimal()

# ggsave(filename = "no_complexes_bar.pdf", 
#        plot = last_plot(), path = "output/F4S1/",
#        width = 4,
#        height = 4)
