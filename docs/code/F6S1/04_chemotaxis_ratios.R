source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/functional/from_this_paper/vMIPII_chemotaxis_expts/vmipii_tcell_chemotaxis_ratios.csv") %>%
  group_by(condition) %>%
  dplyr::summarise(mean = mean(area_hour_4, na.rm=T)) %>%
  ungroup()

cxcl12 <- data %>%
  filter(condition == "cxcl12") %>%
  dplyr::select(mean)
cxcl12 <- cxcl12$mean
vmipii_cxcl12 <- data %>%
  filter(condition == "vmipii_cxcl12") %>%
  dplyr::select(mean)
vmipii_cxcl12 <- vmipii_cxcl12$mean

data <- data %>%
  mutate(fold_change = (cxcl12 - mean)/(cxcl12-vmipii_cxcl12),
         logfold_change = log((cxcl12 - mean)/(cxcl12-vmipii_cxcl12)))
# The numerator compares the reduction of chemotaxis in the presence of a 
# vMIP-II variant to that will CXCL12 (full/maximal chemotaxis). The requation 
# then takes the ratio of the reduction in chemotaxis to the "gain" exhibited 
# by the assay - i.e. the difference between maximal chemotaxis (CXCL12 alone)
# and minimal chemotaxis (CXCL12 PLUS WT vMIP-II). Overall the equation is 
# comparing how much chemotaxis occurs with the vMIP-II mutants versus
# with WT vMIP-II, with most vMIP-II mutants "allowing" a lot of chemotaxis,
# causing a higher numerator and a smaller fraction. When the log is then taken,
# this causes "poor" vMIP-II antyagonist variants to have more negative log
# values
# The closer fold change is to 1.0 (i.e. the higher the numerator), the lower
# the more inhibition there is and the more the vMIP-II mutant resembles WT
# vMIP-II as a "good" CXCR4 antagonist.


data$condition <- factor(data$condition, 
                         levels = c("r7i_cxcl12","k10t_cxcl12","l13f_cxcl12","triple_cxcl12","cxcl12","vmipii_cxcl12"))

data %>%
  filter(!(condition %in% c("cxcl12", "vmipii_cxcl12"))) %>%
  ggplot(aes(condition, logfold_change)) +
  geom_segment(aes(x=condition, xend=condition, y=0, yend=logfold_change)) +
  geom_point(shape = 21, fill = "white", size = 5, stroke = 0.5) +
  ylim(-7,7)+
  theme_minimal()

# ggsave(filename = "vmipii_chemotaxis_ratio.pdf",
#        plot = last_plot(), path = "output/F6/",
#        width = 3,
#        height = 4)
