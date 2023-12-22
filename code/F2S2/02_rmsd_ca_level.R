source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/structure/processed/pairwise_complex_rmsd.csv")
temp <- data %>% dplyr::count(gnccn) %>% dplyr::filter(n >= 60) 
temp <- temp$gnccn
 # count number of pairwise comparisons in which each position is represented;
 # remove CCN found in less than 50% pairwise comparisons (below)
data <- data %>% group_by(gnccn) %>% 
  filter(gnccn %in% c(temp)) %>%
  dplyr::mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  dplyr::select(-RMSD, -file1, -file2) %>% unique() 
rm(temp)

# order positions
order.ck <- as.factor(c("NTc.Cm10","NTc.Cm9","NTc.Cm8","NTc.Cm7","NTc.Cm6","NTc.Cm5","NTc.Cm4","NTc.Cm3","NTc.Cm2","NTc.Cm1", "NTc.Cm0", "CX.1", "CX.2", "CX.3", "CX.4", "CX.5","cxb1.1","cxb1.2","cxb1.5","cxb1.6","cxb1.7","cxb1.8","cxb1.10","cxb1.11","cxb1.14","cxb1.15","cxb1.16","B1.1","B1.2","B1.3","B1.4","B1.5","B1.6","B1.7","b1b2.4","b1b2.6","b1b2.7","b1b2.8", "b1b2.9","b1b2.10","b1b2.12","b1b2.13","b1b2.14","b1b2.16", "b1b2.17", "B2.1","B2.2","B2.3","B2.4","B2.5","B2.6","b2b3.2","b2b3.3","b2b3.4","b2b3.11","b2b3.12","B3.1","B3.2","B3.3","B3.4","b3h.1","b3h.2","b3h.3","b3h.4","H.1","H.2","H.3","H.4","H.5","H.6","H.7","H.8","H.9","H.10", "CT"))
data$gnccn <- factor(data$gnccn, levels = (order.ck))

# plot
data %>%
  filter(!(gnccn %in% c("CT"))) %>%
  ggplot() +
  geom_line(aes(gnccn, mean, group=1)) +
  geom_line(aes(gnccn, mean-sd, group=1)) +
  geom_line(aes(gnccn, mean+sd, group=1)) +
  
  coord_flip() +
  theme_minimal()

# ggsave(filename = "ca_rmsd.pdf", 
#        plot = last_plot(), path = "output/F2S2/",
#        width = 3,
#        height = 5)
