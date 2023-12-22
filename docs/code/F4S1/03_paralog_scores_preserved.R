source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

ck <- read_csv("data/sequence/chemokine/processed/CK_CONSERVATION.csv") %>%
  filter(ccn %in% c("cxb1.1", "cxb1.2", "B3.1"))
ck$ccn <- factor(ck$ccn, levels = c("cxb1.1", "cxb1.2", "B3.1"))

ck %>%
  select(protein, ccn, ortho_cons) %>% 
  unique() %>%
  ggplot(aes(ccn, ortho_cons)) +
  # geom_boxplot() +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2) +
  theme_minimal()

# ggsave(filename = "ck_ortholog_cons.pdf", 
#        plot = last_plot(), path = "output/F4S1/",
#        width = 4,
#        height = 3)

ckr <- read_csv("data/sequence/gpcr/processed/CKR_CONSERVATION.csv") %>%
  filter(gn %in% c("gnNTr.Cm2", "gnNTr.Cm3", "gn7x27")) %>%
  filter(!(protein %in% c("ackr1","ackr2","ackr3","ackr4","ccrl2")))
ckr$gn <- factor(ckr$gn, levels = c("gnNTr.Cm2", "gnNTr.Cm3", "gn7x27"))

ckr %>%
  select(protein, gn, ortho_cons) %>% 
  unique() %>%
  ggplot(aes(gn, ortho_cons)) +
  geom_violin(trim = FALSE) +
  stat_summary(fun.y=median, geom="point", shape=23, size=2) +
  theme_minimal()

# ggsave(filename = "ckr_ortholog_cons.pdf", 
#        plot = last_plot(), path = "output/F4S1/",
#        width = 4,
#        height = 3)

