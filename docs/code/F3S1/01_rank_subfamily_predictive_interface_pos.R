source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# CHEMOKINE
ck <- read_csv("data/sequence/chemokine/processed/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
ck$motif <- factor(ck$motif, levels = ck$motif[order(ck$mean)])
  
# subset by interface in chemokine-GPCR complexes
interface <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file != "6meo") %>%
  dplyr::select(source_gnccn) %>% unique()
interface <- unique(interface$source_gnccn)
  
# plot top
ck %>%
  filter(motif %in% interface) %>%
  filter(mean >= .75) %>%
  ggplot(aes(motif, mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
  coord_flip() +
  ylim(.49,1) +
  theme_minimal()
  
# how many
temp <- ck  %>%
  filter(mean >= .75) %>% unique() %>%
  filter(motif %in% interface) %>% unique()
paste0("There are ", nrow(temp), " chemokine subfamily-predictive positions at the chemokine-GPCR interface")  

# ggsave(filename = "chemokine_top_subfamily.pdf", 
#        plot = last_plot(), path = "output/F3S1/",
#        width = 3,
#        height = 8)
  
# RECEPTOR
ckr <- read_csv("data/sequence/gpcr/processed/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>% 
  separate(motif, into = c("temp", "motif"), sep = "gn") %>% dplyr::select(-temp)
ckr$motif <- factor(ckr$motif, levels = ckr$motif[order(ckr$mean)])

# subset by interface in chemokine-GPCR complexes
interface <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file != "6meo") %>%
  dplyr::select(target_gnccn) %>% unique()
interface <- unique(interface$target_gnccn)
  
# plot top
ckr %>%
  filter(motif %in% interface) %>%
  #top_n(20, mean) %>%
  filter(mean >= .75) %>%
  ggplot(aes(motif, mean)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(0.05)) +
  geom_point(shape = 21, colour = "black", fill = "white", size = 4, stroke = 0.5) +
  coord_flip() +
  ylim(.49,1) +
  theme_minimal()
  
# count how many
temp <- ckr %>%
  filter(mean >= .75) %>%
  filter(motif %in% interface) %>% unique()
paste0("There are ", nrow(temp), " chemokine subfamily-predictive positions at the chemokine-GPCR interface")  

# ggsave(filename = "receptor_top_subfamily.pdf", 
#        plot = last_plot(), path = "output/F3S1/",
#        width = 3,
#        height = 8)
  