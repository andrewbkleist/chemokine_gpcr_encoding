source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
# (1) CHEMOKINE ----------------------------------------------------------------
# import rin
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo") %>%
  dplyr::select(source_gnccn, dom1, file)
  
# add labels
rin <- rin %>% mutate(str_ck = case_when(
  dom1 == "NTc" ~ "NTc",
  dom1 != "NTc" ~ "core_c"
))
  
# summary stats
contacts <- rin %>% dplyr::count(source_gnccn, file, str_ck)

# summary stats and p-value
# see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
sum <- wilcox.test(n ~ str_ck, data = contacts, exact = FALSE)
sum$p.value
  
# plot
contacts %>%
  ggplot(aes(str_ck, n)) +
  geom_boxplot() +
  ylim(0,7) +
  theme_minimal() 

# ggsave(filename = "ck_str_unstr_contacts_per_residue_box.pdf",
#        plot = last_plot(), path = "output/F4S2/",
#        width = 3,
#        height = 7)

# (2) GPCR ---------------------------------------------------------------------  
# import rin
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(file != "6meo") %>%
  dplyr::select(target_gnccn, dom2, file)
  
# add labels
rin <- rin %>% mutate(str_ckr = case_when(
  dom2 == "NTr" ~ "NTr",
  dom2 != "NTr" ~ "core_r"
))
  
# summary stats and p-value
# see http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r
contacts <- rin %>% dplyr::count(target_gnccn, file, str_ckr)

ckr.nt_core <- filter(contacts, str_ckr == "NTr" | str_ckr == "core_r" )
ckr.nt_core <- wilcox.test(n ~ str_ckr, data = ckr.nt_core, exact = FALSE)
ckr.nt_core$p.value # significant
  
# plot
contacts %>%
  ggplot(aes(str_ckr, n)) +
  geom_boxplot() +
  ylim(0,7) +
  theme_minimal() 

# ggsave(filename = "ckr_str_unstr_contacts_per_residue_box.pdf",
#        plot = last_plot(), path = "output/F4S2/",
#        width = 3,
#        height = 7)
