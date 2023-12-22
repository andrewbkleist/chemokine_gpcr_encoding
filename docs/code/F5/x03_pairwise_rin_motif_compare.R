source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import ck motifs
# motif.ck.nt <- read_csv("data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CK_NT.csv") %>%
#   separate(col = pair, into = c("ck1", "ck2"), sep = "_") %>%
#   dplyr::count(ck1, ck2, type) %>%
#   group_by(ck1, ck2) %>%
#   mutate(sum_ck = sum(n)) %>%
#   ungroup() %>%
#   mutate(ck_pid = n/sum_ck) %>%
#   filter(type == "a_and_b") %>%
#   dplyr::select(-n, -sum_ck, -type)
motif.ck.nt <- read_csv("data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CK_NT.csv") %>%
  separate(col = pair, into = c("ck1", "ck2"), sep = "_") %>%
  dplyr::count(ck1, ck2, type) %>%
  pivot_wider(names_from = type, values_from = n)
motif.ck.nt[is.na(motif.ck.nt)] <- 0
motif.ck.nt <- motif.ck.nt %>%
  pivot_longer(cols = c("a_and_b","a_not_b", "b_not_a"), names_to = "type") %>%
  group_by(ck1, ck2) %>%
  mutate(sum_ck = sum(value)) %>%
  ungroup() %>%
  mutate(ck_pid = value/sum_ck) %>%
  filter(type == "a_and_b") %>%
  dplyr::select(-value, -sum_ck, -type)
  
# import ckr nt motifs
# motif.ckr.nt <- read_csv("data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CKR_NT.csv") %>%
#   separate(col = pair, into = c("ckr1", "ckr2"), sep = "_") %>%
#   dplyr::count(ckr1, ckr2, type) %>%
#   group_by(ckr1, ckr2) %>%
#   mutate(sum_ckr = sum(n)) %>%
#   ungroup() %>%
#   mutate(ck_pid = n/sum_ckr) %>%
#   filter(type == "a_and_b") %>%
#   dplyr::select(-n, -sum_ckr, -type)
motif.ckr.nt <- read_csv("data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CKR_NT.csv") %>%
  separate(col = pair, into = c("ckr1", "ckr2"), sep = "_") %>%
  dplyr::count(ckr1, ckr2, type) %>%
  pivot_wider(names_from = type, values_from = n)
motif.ckr.nt[is.na(motif.ckr.nt)] <- 0
motif.ckr.nt <- motif.ckr.nt %>%
  pivot_longer(cols = c("a_and_b","a_not_b", "b_not_a"), names_to = "type") %>%
  group_by(ckr1, ckr2) %>%
  mutate(sum_ckr = sum(value)) %>%
  ungroup() %>%
  mutate(ckr_pid_nt = value/sum_ckr) %>%
  filter(type == "a_and_b") %>%
  dplyr::select(-value, -sum_ckr, -type)

# import ckr ecl2 motifs
# motif.ckr.ecl2 <- read_csv("data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CKR_ECL2.csv") %>%
#   separate(col = pair, into = c("ckr1", "ckr2"), sep = "_") %>%
#   dplyr::count(ckr1, ckr2, type) %>%
#   group_by(ckr1, ckr2) %>%
#   mutate(sum_ckr = sum(n)) %>%
#   ungroup() %>%
#   mutate(ck_pid = n/sum_ckr) %>%
#   filter(type == "a_and_b") %>%
#   dplyr::select(-n, -sum_ckr, -type)
motif.ckr.ecl2 <- read_csv("data/motif/processed/PAIRWISE_MOTIF_COMPARISONS_CKR_ECL2.csv") %>%
  separate(col = pair, into = c("ckr1", "ckr2"), sep = "_") %>%
  dplyr::count(ckr1, ckr2, type) %>%
  pivot_wider(names_from = type, values_from = n)
motif.ckr.ecl2[is.na(motif.ckr.ecl2)] <- 0
motif.ckr.ecl2 <- motif.ckr.ecl2 %>%
  pivot_longer(cols = c("a_and_b","a_not_b", "b_not_a"), names_to = "type") %>%
  group_by(ckr1, ckr2) %>%
  mutate(sum_ckr = sum(value)) %>%
  ungroup() %>%
  mutate(ckr_pid_ecl2 = value/sum_ckr) %>%
  filter(type == "a_and_b") %>%
  dplyr::select(-value, -sum_ckr, -type)

# import structure
rin <- read_csv("data/structure/processed/RIN_all_by_all_pid.csv")
colnames(rin) <- c("pdb1", "pdb2", "rin_pid", "ck1", "ckr1", "ck2", "ckr2")
rin <- left_join(rin, motif.ck.nt)
rin <- left_join(rin, motif.ckr.nt)
rin <- left_join(rin, motif.ckr.ecl2)
rin$motif_pid_ave <- rowMeans(subset(rin, select = c(ck_pid, ckr_pid_nt, ckr_pid_ecl2)), na.rm = FALSE)

# remove complexes in which non-native or non-human sequences are used 
# native CCL5-CCR5 (7f1r) excluded since it is incomplete from structural POV
rin <- rin %>%
  filter(!(pdb1 %in% c("5uiw","7o7f","4rws", "4xt1","5wb2", "7f1r"))) %>%
  filter(!(pdb2 %in% c("5uiw","7o7f","4rws", "4xt1","5wb2","7f1r"))) 

# unite columns
rin <- rin %>% unite(ck1_ck2, c(ck1, ck2))
rin <- rin %>% unite(ckr1_ckr2, c(ckr1, ckr2))

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
rin$type1 <- lookup$type[match(unlist(rin$pdb1), lookup$pdb)]
rin$type2 <- lookup$type[match(unlist(rin$pdb2), lookup$pdb)]
rin <- rin %>% unite(type12, c(type1, type2))
rin <- rin %>% 
  mutate(type12summ = case_when(
    type12 %in% c("ccl5_rltd_ccl5_rltd", "cxcl8_rltd_cxcl8_rltd") ~ "intra_network_same_sub",
    type12 %in% c("other_cxc_cxcl8_rltd", "other_cc_ccl5_rltd") ~ "inter_network_same_sub",
    type12 %in% c("cxc_ccl5_rltd", "other_other_cc", "cxcl8_rltd_ccl5_rltd", "other_ccl5_rltd", "cxcl8_rltd_other_cc", "other_other", "other_cxc_ccl5_rltd", "other_cxc_other_cc", "other_cxcl8_rltd", "other_other_cxc") ~ "diff_sub",
    type12 %in% c("viral_ccl5_rltd", "viral_other_cc", "viral_cxc", "viral_other", "viral_viral", "viral_cxcl8_rltd", "viral_other_cxc") ~ "viral"
  ))

# plot
rin %>%
  ggplot(aes(ckr_pid_ecl2, rin_pid, fill = type12summ)) +
  geom_jitter(shape = 21, colour = "black",  size = 9, stroke = 0.5) +
  scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  xlim(0,0.05)+
  
  theme_minimal()

# ggsave(filename = "pairwise_rin_seq_compare.pdf",
#        plot = last_plot(), path = "output/F5/",
#        width = 6,
#        height = 4)

p <- rin %>%
  ggplot(aes(motif_pid_ave, rin_pid, fill = ck1_ck2)) +
  geom_jitter(shape = 21, colour = "black",  size = 9, stroke = 0.5) +
  # scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  # xlim(0,0.05)+
  theme_minimal()
p
library(plotly)
ggplotly(p)


