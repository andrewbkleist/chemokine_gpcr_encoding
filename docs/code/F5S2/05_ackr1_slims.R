source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# (1) ACKR1 FINGERPTINTS FYA VS FYB --------------------------------------------
# import, reformat
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data <- data %>%
  filter(motif %in% c("YD", "DYD", "YDA", "GDYD", "DYDA", "YDAN",
                      "YG", "DYG", "YGA", "GDYG", "DYGA", "YGAN")) %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  mutate(fya_fyb = case_when(
    motif %in% c("YD", "DYD", "YDA", "GDYD", "DYDA", "YDAN") ~ "d42",
    motif %in% c("YG", "DYG", "YGA", "GDYG", "DYGA", "YGAN") ~ "g42"
  ))

# order receptors
order.ckr <- as.factor(tolower(c("ccr1","ccr2", "ccr3", "ccr4",
                                 "ccr5","ccr6","ccr7","ccr8","ccr9",
                                 "ccr10","cxcr1","cxcr2","cxcr3",
                                 "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                 "ackr1","ackr2","ackr3","ackr4",
                                 "ccrl2")))
levels(data$protein)
data$protein <- factor(data$protein, levels = order.ckr)

# select MOTIF "SLICE"
order.mot <- as.factor(unique(c("YD", "DYD", "YDA", "GDYD", "DYDA", "YDAN",
                                "YG", "DYG", "YGA", "GDYG", "DYGA", "YGAN")))
data$motif <- factor(data$motif, levels = rev(order.mot))

data %>%
  
  ggplot() + 
  geom_tile(aes(protein, motif), fill = "mediumpurple3") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))

# ggsave(filename = "ackr1_slims_fingerprint.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 5,
#        height = 3)



# (2) ACKR1 ORTHOLOG CONS BOXPLOTS ---------------------------------------------

data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data <- data %>%
  filter(motif %in% c("YD", "DYD", "YDA", "GDYD", "DYDA", "YDAN",
                      "YG", "DYG", "YGA", "GDYG", "DYGA", "YGAN")) %>%  
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  mutate(fya_fyb = case_when(
    motif %in% c("YD", "DYD", "YDA", "GDYD", "DYDA", "YDAN") ~ "d42",
    motif %in% c("YG", "DYG", "YGA", "GDYG", "DYGA", "YGAN") ~ "g42"
  ))
data %>%
  filter(protein != "ackr1") %>%
  ggplot(aes(fya_fyb, pct_ortho)) +
  geom_boxplot() +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 8, stroke = 0.5) +
  # stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
  #              geom = "crossbar", width = 0.5) +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  coord_flip() +
  theme_minimal()

# ggsave(filename = "ackr1_slims_box.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 7,
#        height = 3)

# stats
# http://www.sthda.com/english/wiki/normality-test-in-r#google_vignette
shapiro.test(data$pct_ortho)
# p-value = 0.04561 (less than 0.05) therefore must use non-parametric
# http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
test <- data %>% filter(protein != "ackr1")
test <- wilcox.test(pct_ortho ~ fya_fyb, data = data,
                   exact = FALSE)
test




# OLD, NOT USED
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter(!(motif == "DYG" & protein == "ackr1")) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.55) +
  theme_minimal()

# ggsave(filename = "ackr1_dyg_ortho_para.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 4.5,
#        height = 4)

data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter((motif == "DYG" & protein == "ackr1") |  (motif == "DYD" & protein == "cxcr4") ) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  ylim(0,1) +
  xlim(0,0.55) +
  theme_minimal()

# ggsave(filename = "ackr1_dyg2_ortho_para.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 4.5,
#        height = 4)



################

data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  filter( (motif == "DYG" & protein == "ackr1") |  motif == "DYD" ) %>%
  # filter(!(motif != "GMP" & protein == "cxcr1")) %>%
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.01, height = 0.01, shape = 21, colour = "white", size = 8, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  ylim(0,1) +
  # xlim(0,0.55) +
  theme_minimal()

# ggsave(filename = "ckr_nt_motif_ortho_para_ackr1.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 5,
#        height = 4)


#---
# EXPLORATION 20231209

# DOTPLOT
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  # filter(motif %in% c("DY", "YD", "DYD", "EY","YE", "DYG")) %>%
  filter(motif %in% c("YD", "YG")) %>%
  
  ggplot(aes(pct_super, pct_ortho)) +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  # scale_fill_gradient2(low = "white", high = "mediumpurple3") +
  scale_fill_manual(values=c("gray70", "mediumpurple3")) + 
  # ylim(0,1) +
  # xlim(0,0.55) +
  theme_minimal()

# COMPARISON 
data <- read_csv("data/motif/processed/CKR_MOTIF_FREQUENCY_NTERM.csv")
data <- data %>%
  filter(motif %in% c("YD", "YG", "DYD", "DYG", "YDA","YGA", "DYGA","DYDA", "YGAN","YDAN")) %>%
  mutate(slim = case_when(
    pct_ortho >= 0.5 ~ "yes",
    pct_ortho < 0.5 ~ "no"
  )) %>%
  mutate(fya_fyb = case_when(
    motif %in% c("YD", "DYD", "YDA", "DYDA", "YGAN") ~ "d42",
    motif %in% c("YG", "DYG", "YGA", "DYGA", "YDAN") ~ "g42"
  ))
data %>%
  filter(protein != "ackr1") %>%
  ggplot(aes(fya_fyb, pct_ortho)) +
  geom_boxplot() +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  theme_minimal()

data %>%
  ggplot(aes(fya_fyb, count_super)) +
  geom_boxplot() +
  geom_jitter(aes(fill = slim), width = 0.05, height = 0.05, shape = 21, colour = "white", size = 6, stroke = 0.5) +
  theme_minimal()

# ADD TILE PLOT/FINGERPRINT

# (1) SULFOTYROSINE MOTIFS


