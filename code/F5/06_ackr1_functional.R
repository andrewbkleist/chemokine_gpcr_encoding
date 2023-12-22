source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/functional/from_this_paper/ACKR1_expts/ackr1_ec50s.csv") %>%
  pivot_wider(names_from = ackr1, values_from = ec50) %>%
  mutate(log2ratio = log2(wt/g42d))

data$chemokine <- factor(data$chemokine, levels = c("ccl2","ccl7","cxcl1","cxcl8","cxcl11","cxcl12"))
# log2ratio
data %>%
  filter(chemokine != ("cxcl12")) %>%
  ggplot(aes(chemokine, log2ratio)) +
  geom_segment(aes(x=chemokine, xend=chemokine, y=0, yend=log2ratio)) +
  geom_point(shape = 21, fill = "white", size = 5, stroke = 0.5) +
  ylim(-1,2)+
  theme_minimal()

# ggsave(filename = "ackr1_log2ec50.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 3,
#        height = 3)


