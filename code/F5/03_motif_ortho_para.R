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
  filter(motif != "KPV") %>%
  filter(!(motif != "ELR" & protein == "cxcl12")) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.25) +
  theme_minimal()

# ggsave(filename = "ck_motif_ortho_para.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 4.5,
#        height = 4)

data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter((motif == "KPV" & protein == "cxcl12") | (motif == "ELR" & protein == "cxcl8")) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.25) +
  theme_minimal()

# ggsave(filename = "ck_motif_ortho_para_elr_kpv.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 4.5,
#        height = 4)

# (2) GPCR NT ------------------------------------------------------------------

data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter(!(motif != "DYG" & protein == "ccr1")) %>%
  filter(!(motif != "GMP" & protein == "cxcr1")) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.55) +
  theme_minimal()

# ggsave(filename = "ckr_nt_motif_ortho_para.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 4.5,
#        height = 4)

data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter((motif == "DYG" & protein == "ccr1") | (motif == "GMP" & protein == "cxcr1")) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.55) +
  theme_minimal()

# ggsave(filename = "ckr_nt_motif_ortho_para_dyg_gmp.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 4.5,
#        height = 4)
