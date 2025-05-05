#

# packages, working directory
library(tidyverse)
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")

################################################################################


# mer3
data <- read_csv("F5/020117_CCR10_CCR3_Nterm/integrated_ec50_conservation.csv")  

data$motif <- factor(data$motif, levels = c("SEA","EAI","AIL","ILP","LPI","PIA","IAS","ASS","SSC"))
data$receptor <- factor(data$receptor, levels = c("ccr3", "ccr10"))

data %>%
  filter(motif_len == "mer3") %>%
  filter(motif %in% c("SEA", "EAI", "AIL", "ILP", "LPI")) %>%
  ggplot(aes(motif,logec50)) +
  geom_path(group=1) +     
  geom_errorbar(aes(ymin=logec50-logec50error, ymax=logec50+logec50error), width=.3,
                position=position_dodge(0.05)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
  theme_minimal() +
  coord_flip() +
  facet_grid(.  ~ receptor)

ggsave(filename = "ccl28_functional.pdf", 
       plot = last_plot(), path = "F5/output/",
       width = 3.5,
       height = 4)

data %>%
  filter(motif_len == "mer3") %>%
  filter(motif %in% c("SEA", "EAI", "AIL", "ILP", "LPI")) %>%
  ggplot(aes(motif,ortholog_cons)) +
  geom_path(group=1) +     
  geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
  theme_minimal() +
  ylim(0,1) +
  coord_flip() +
  facet_grid(. ~ receptor)

ggsave(filename = "ccl28_conservation.pdf", 
       plot = last_plot(), path = "F5/output/",
       width = 3.5,
       height = 4)



# mer2
data <- read_csv("F5/020117_CCR10_CCR3_Nterm/integrated_ec50_conservation.csv")  

data$motif <- factor(data$motif, levels = c("SE","EA","AI","IL","LP","PI","IA","AS","SS"))
data$receptor <- factor(data$receptor, levels = c("ccr3", "ccr10"))

data %>%
  filter(motif_len == "mer2") %>%
  ggplot(aes(motif,logec50)) +
  geom_path(group=1) +     
  geom_errorbar(aes(ymin=logec50-logec50error, ymax=logec50+logec50error), width=.3,
                position=position_dodge(0.05)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
  theme_minimal() +
  facet_grid(receptor ~ .)

# ggsave(filename = "ccl28_functional.pdf", 
#        plot = last_plot(), path = "F5/output/",
#        width = 6,
#        height = 6)

data %>%
  filter(motif_len == "mer2") %>%
  ggplot(aes(motif,ortholog_cons)) +
  geom_path(group=1) +     
  geom_point(shape = 21, colour = "black", fill = "white", size = 3, stroke = 0.5) +
  theme_minimal() +
  facet_grid(receptor ~ .)
