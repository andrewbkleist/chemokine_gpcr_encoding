source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CHEMOKINE NT -------------------------------------------------------------
data <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif, count_super) %>%
  unique()
data$motif <- factor(data$motif, levels = data$motif[order(data$count_super)])

# plot top
data %>%
  top_n(20, count_super) %>%
  ggplot(aes(count_super, motif)) +
  # geom_bar(stat = "identity") +
  geom_point(shape = 21, colour = "black", fill = "mediumpurple3", size = 5, stroke = 0.5) +
  # coord_flip() +
  xlim(0,11) +
  theme_minimal()

ggsave(filename = "top_slims_ck_nt.pdf",
       plot = last_plot(), path = "output/F5S1/",
       width = 3,
       height = 8)

# (2) GPCR NT ------------------------------------------------------------------
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif, count_super) %>%
  unique()
data$motif <- factor(data$motif, levels = data$motif[order(data$count_super)])

# plot top
data %>%
  top_n(20, count_super) %>%
  ggplot(aes(count_super, motif)) +
  # geom_bar(stat = "identity") +
  geom_point(shape = 21, colour = "black", fill = "mediumpurple3", size = 5, stroke = 0.5) +
  # coord_flip() +
  xlim(0,15) +
  theme_minimal()

ggsave(filename = "top_slims_ckr_nt.pdf",
       plot = last_plot(), path = "output/F5S1/",
       width = 3,
       height = 8)

# (2) GPCR ECL2 ----------------------------------------------------------------
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_ECL2.csv") %>%
  filter(pct_ortho>= 0.5) %>%
  dplyr::select(motif, count_super) %>%
  unique()
data$motif <- factor(data$motif, levels = data$motif[order(data$count_super)])

# plot top
data %>%
  top_n(20, count_super) %>%
  ggplot(aes(count_super, motif)) +
  # geom_bar(stat = "identity") +
  geom_point(shape = 21, colour = "black", fill = "mediumpurple3", size = 5, stroke = 0.5) +
  # coord_flip() +
  xlim(0,7) +
  theme_minimal()

ggsave(filename = "top_slims_ckr_ecl2.pdf",
       plot = last_plot(), path = "output/F5S1/",
       width = 3,
       height = 8)

