source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CCL28 -----------------------------------------------------------------------
data <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif, count_super) %>% 
  unique() %>%
  dplyr::count(count_super)
total <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif) %>% 
  unique()
data <- data %>% mutate(pct = n/nrow(total))
  
# graph
hsize <- 1
data %>%
  ggplot(aes(x = 1, y= pct)) +
  geom_col() +
  coord_polar("y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()

# ggsave(filename = "pie_pct_uniqueness_ck_nt.pdf",
#        plot = last_plot(), path = "output/F5S1/",
#        width = 3,
#        height = 3)

# (2) RECEPTOR NT --------------------------------------------------------------
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif, count_super) %>% 
  unique() %>%
  dplyr::count(count_super)
total <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif) %>% 
  unique()
data <- data %>% mutate(pct = n/nrow(total))

# graph
hsize <- 1
data %>%
  ggplot(aes(x = 1, y= pct)) +
  geom_col() +
  coord_polar("y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()

ggsave(filename = "pie_pct_uniqueness_ckr_nt.pdf",
       plot = last_plot(), path = "output/F5S1/",
       width = 3,
       height = 3)


# (3) RECEPTOR ECL2 ------------------------------------------------------------
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif, count_super) %>% 
  unique() %>%
  dplyr::count(count_super)
total <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif) %>% 
  unique()
data <- data %>% mutate(pct = n/nrow(total))

# graph
hsize <- 1
data %>%
  ggplot(aes(x = 1, y= pct)) +
  geom_col() +
  coord_polar("y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()

ggsave(filename = "pie_pct_uniqueness_ckr_ecl2.pdf",
       plot = last_plot(), path = "output/F5S1/",
       width = 3,
       height = 3)

  