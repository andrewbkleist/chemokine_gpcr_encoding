source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) CHEMOKINE NT -------------------------------------------------------------
data <- read_csv("data/motif/processed/CK_MOTIF_FREQUENCY.csv")
data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter(motif %in% c("SEA", "EAI", "AIL", "ILP", "LPI", "PIA", "IAS", "ASS")) %>%
  # filter(motif %in% c("SE", "EA", "AI", "IL", "LP", "PI")) %>%
  filter(protein == "ccl28") %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_point(aes(fill = slim), width = 0.005, height = 0.005, shape = 21, colour = "white", size = 9, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("grey70", "mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.25) +
  theme_minimal()
temp <- data %>%   filter(motif %in% c("SEA", "EAI", "AIL", "ILP", "LPI", "PIA", "IAS", "ASS")) %>%
  filter(protein == "ccl28")

# ggsave(filename = "ckr_nt_motif_ortho_para_ccl28.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 5,
#        height = 4)
