source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# read csv
df <- read_csv("data/structure/processed/RIN_all_by_all_pid.csv")

# unite columns
df <- df %>% unite(ck1_ck2, c(ck1, ck2))
df <- df %>% unite(ckr1_ckr2, c(ckr1, ckr2))
rm(lookup)

# map pdb "classes"
ccl5_rltd <- data.frame(pdb = c("5uiw", "7o7f", "7f1r", "zheng", "7vl9", "7xa3", "7f1t"), type = "ccl5_rltd")
other_cc <- data.frame(pdb = c("6wwz"), type = "other_cc")
cxcl8_cxc <- data.frame(pdb = c("6lfo", "8ic0"), type = "cxcl8_rltd")
other_cxc <- data.frame(pdb = c("ngo"), type = "other_cxc")
other <- data.frame(pdb = c("7sk3", "7xbx"), type = "other")
viral <- data.frame(pdb = c("4rws", "4xt1", "5wb2"), type = "viral")
lookup <- rbind(ccl5_rltd,other_cc, cxcl8_cxc, other_cxc, other, viral)
rm(ccl5_rltd,other_cc, cxcl8_cxc, other_cxc, other, viral)

# define pairwise relationships
df$type1 <- lookup$type[match(unlist(df$a), lookup$pdb)]
df$type2 <- lookup$type[match(unlist(df$b), lookup$pdb)]
df <- df %>% unite(type12, c(type1, type2))
df <- df %>% 
  mutate(type12summ = case_when(
    type12 %in% c("ccl5_rltd_ccl5_rltd", "cxcl8_rltd_cxcl8_rltd") ~ "intra_network_same_sub",
    type12 %in% c("other_cxc_cxcl8_rltd", "other_cc_ccl5_rltd") ~ "inter_network_same_sub",
    type12 %in% c("cxc_ccl5_rltd", "other_other_cc", "cxcl8_rltd_ccl5_rltd", "other_ccl5_rltd", "cxcl8_rltd_other_cc", "other_other", "other_cxc_ccl5_rltd", "other_cxc_other_cc", "other_cxcl8_rltd", "other_other_cxc") ~ "diff_sub",
    type12 %in% c("viral_ccl5_rltd", "viral_other_cc", "viral_cxc", "viral_other", "viral_viral", "viral_cxcl8_rltd", "viral_other_cxc") ~ "viral"
  ))

# define median pairwise
median(df$n)
paste0("The mean number of shared contacts among any two complexes is ",
       round(median(df$n)*100, 2), "%")

# plot
df$arb <- c(1)
df %>%
  filter(type12summ != "viral") %>%
  ggplot(aes(arb, n, fill = type12summ)) +
  geom_jitter(shape = 21, colour = "black",  size = 6, stroke = 0.5) +
  ylim(-0.05, 0.52) +
  geom_hline(yintercept = c(0.1, 0.2)) +
  scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  theme_minimal()

# ggsave(filename = "rin_pct_compare_dot.pdf", 
#        plot = last_plot(), path = "output/F4/",
#        width = 3.5,
#        height = 6)

# df %>%
#   filter(type12summ != "viral") %>%
#   ggplot(aes(x = n, y = type12summ, fill = type12summ)) +
#   geom_density_ridges() +
#   scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
#   theme_minimal() 
