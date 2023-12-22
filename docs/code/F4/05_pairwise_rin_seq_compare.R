# Imports chemokine and GPCR pairwise sequence comparisons and calculates 
# sequence PID based on alignments. In order to focus on sequence identity at 
# positions that could influence selectivity, only positions that make contacts
# in any of the 16 complexes are considered. Note that all interface positions
# are included for all chemokines/GPCRs, e.g. even if a particular position
# does not make a contact in an existing structure, that position is still
# considered for the sequence PID calculation as a position that could
# theoretically make a contact with another chemokine or GPCR

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import ck sequences
interface.ck <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("source_gnccn") %>% unique() %>%
  filter(!is.na(source_gnccn))
interface.ck <- interface.ck$source_gnccn
seq.ck <- read_csv("data/sequence/chemokine/processed/CK_ALL_BY_ALL_PAIRWISE_IDENTITY_BY_POS.csv") %>%
  separate(col = pair, into = c("ck1", "ck2"), sep = "_") %>%
  mutate(ck_pid = identity) %>%
  filter(!(resid_a == "-" & resid_b == "-")) %>%
  dplyr::select(-resid_a, -resid_b, -identity) %>%
  filter(source_gnccn %in% c(interface.ck)) %>%
  dplyr::count(ck1, ck2, ck_pid) %>%
  group_by(ck1, ck2) %>%
  mutate(sum_ck = sum(n)) %>%
  ungroup() %>%
  filter(ck_pid == 1) %>%
  mutate(ck_pid = n/sum_ck) %>%
  dplyr::select(-n, -sum_ck)
rm(interface.ck)

# import ckr sequences
interface.ckr <- read_csv("data/structure/processed/RIN_residue.csv") %>% 
  filter(Chain1 != Chain2) %>%
  dplyr::select("target_gnccn") %>% unique() %>%
  filter(!is.na(target_gnccn))
interface.ckr <- interface.ckr$target_gnccn
seq.ckr <- read_csv("data/sequence/gpcr/processed/CKR_ALL_BY_ALL_PAIRWISE_IDENTITY_BY_POS.csv") %>%
  separate(col = pair, into = c("ckr1", "ckr2"), sep = "_") %>%
  mutate(ckr_pid = identity) %>%
  filter(!(resid_a == "-" & resid_b == "-")) %>%
  dplyr::select(-resid_a, -resid_b, -identity) %>%
  filter(target_gnccn %in% c(interface.ckr)) %>%
  dplyr::count(ckr1, ckr2, ckr_pid) %>%
  group_by(ckr1, ckr2) %>%
  mutate(sum_ckr = sum(n)) %>%
  ungroup() %>%
  filter(ckr_pid == 1) %>%
  mutate(ckr_pid = n/sum_ckr) %>%
  dplyr::select(-n, -sum_ckr)
rm(interface.ckr)

# import structure
rin <- read_csv("data/structure/processed/RIN_all_by_all_pid.csv")
colnames(rin) <- c("pdb1", "pdb2", "rin_pid", "ck1", "ckr1", "ck2", "ckr2")
rin <- left_join(rin, seq.ck)
rin <- left_join(rin, seq.ckr)

# remove complexes in which non-native or non-human sequences are used 
# native CCL5-CCR5 (7f1r) excluded since it is incomplete from structural POV
rin <- rin %>%
  filter(!(pdb1 %in% c("5uiw","7o7f","4rws", "4xt1","5wb2", "7f1r"))) %>%
  filter(!(pdb2 %in% c("5uiw","7o7f","4rws", "4xt1","5wb2","7f1r"))) 

# calculate average sequence PID
rin$seq_pid <- (rin$ck_pid + rin$ckr_pid)/2

#--
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
  ggplot(aes(seq_pid, rin_pid, fill = type12summ)) +
  geom_point(shape = 21, colour = "black",  size = 9, stroke = 0.5) +
  scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  theme_minimal()

ggsave(filename = "pairwise_rin_seq_compare.pdf",
       plot = last_plot(), path = "output/F4/",
       width = 6,
       height = 4)

p <- rin %>%
  ggplot(aes(seq_pid, rin_pid, fill = ck1_ck2)) +
  geom_point(shape = 21, colour = "black",  size = 6, stroke = 0.5) +
  # scale_fill_manual(values=c("grey30",  "steelblue3",  "mediumpurple3")) + # "slateblue"/ "dogerblue4" also works
  theme_minimal()
p
library(plotly)
ggplotly(p)


