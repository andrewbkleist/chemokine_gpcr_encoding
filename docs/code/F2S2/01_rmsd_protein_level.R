source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

a <- read_csv("data/structure/processed/pairwise_chemokine_rmsd.csv")
b <- read_csv("data/structure/processed/pairwise_receptor_rmsd.csv")
ab <- read_csv("data/structure/processed/pairwise_complex_rmsd.csv")

a <- a %>% group_by(file1, file2) %>% 
  dplyr::mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  dplyr::select(-gnccn, -RMSD) %>% unique()
a$type <- c("ck")

b <- b %>% group_by(file1, file2) %>% 
  dplyr::mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  dplyr::select(-gnccn, -RMSD) %>% unique()
b$type <- c("ckr")

ab <- ab %>% group_by(file1, file2) %>% 
  dplyr::mutate(mean = mean(RMSD), sd = sd(RMSD)) %>% ungroup() %>% 
  dplyr::select(-gnccn, -RMSD) %>% unique()
ab$type <- c("complex")

data <- rbind(a, b, ab)
rm(a,b,ab)

data %>%
   # filter(file1 != "7sk3") %>% 
   # filter(file2 != "7sk3") %>%
  ggplot(aes(type, mean)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2) +
  theme_minimal()

# ggsave(filename = "fullcomplex_rmsd_violin.pdf",
#        plot = last_plot(), path = "output/F2S2/",
#        width = 3,
#        height = 5)

stats <- data %>% 
  group_by(type) %>%
  dplyr::summarise(median = median(mean)) %>%
  ungroup()
stats
