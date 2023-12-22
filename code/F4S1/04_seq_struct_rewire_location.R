source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# chemokine
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  # filter(!(file %in% c("6meo", "7sk3"))) %>%
  filter(!(file %in% c("5uiw","7o7f","4rws", "4xt1","5wb2", "7f1r", "6meo"))) 
temp <- rin %>%
  dplyr::count(source_gnccn, target_gnccn)
rin <- left_join(rin, temp)
rm(temp)
rin <- rin %>%
  filter(all_para_ck < 0.5) %>%
  filter(all_non_ackr_para_ckr < 0.5) %>%
  filter(n <5) %>% # bottom 25%
  dplyr::count(dom1) %>%
  mutate(total = sum(n)) %>%
  mutate(pct_identity = n / total)
rin$dom1 <- factor(rin$dom1, levels = rin$dom1[order(rin$pct_identity)])
rin %>%
  ggplot(aes(dom1, pct_identity)) +
  geom_bar(stat = "identity") +
  # geom_line(aes(x = dom, y = 1-pct_identity, group = 1)) +
  # geom_point(shape = 21, colour = "white", fill = "black", size = 5, stroke = 2) +
  ylim(0,0.5)+
  coord_flip() +
  theme_minimal()

# ggsave(filename = "rewiring_ck.pdf",
#        plot = last_plot(), path = "output/F4S1/",
#        width = 4,
#        height = 3)

# receptor
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  # filter(!(file %in% c("6meo", "7sk3"))) %>%
  filter(!(file %in% c("5uiw","7o7f","4rws", "4xt1","5wb2", "7f1r", "6meo"))) 
temp <- rin %>%
  dplyr::count(source_gnccn, target_gnccn)
rin <- left_join(rin, temp)
rm(temp)
rin <- rin %>%
  filter(all_para_ck < 0.5) %>%
  filter(all_non_ackr_para_ckr < 0.5) %>%
  filter(n <5) %>% # bottom 25%
  dplyr::count(dom2) %>%
  mutate(total = sum(n)) %>%
  mutate(pct_identity = n / total)
rin$dom2 <- factor(rin$dom2, levels = rin$dom2[order(rin$pct_identity)])
rin %>%
  ggplot(aes(dom2, pct_identity)) +
  geom_bar(stat = "identity") +
  # geom_line(aes(x = dom, y = 1-pct_identity, group = 1)) +
  # geom_point(shape = 21, colour = "white", fill = "black", size = 5, stroke = 2) +
  ylim(0,0.5)+
  coord_flip() +
  theme_minimal()

# ggsave(filename = "rewiring_ckr.pdf",
#        plot = last_plot(), path = "output/F4S1/",
#        width = 4,
#        height = 3)

# joint
rin <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  # filter(!(file %in% c("6meo", "7sk3"))) %>%
  filter(!(file %in% c("5uiw","7o7f","4rws", "4xt1","5wb2", "7f1r", "6meo"))) 
temp <- rin %>%
  dplyr::count(source_gnccn, target_gnccn)
rin <- left_join(rin, temp)
rm(temp)
rin <- rin %>%
  filter(all_non_ackr_para_ckr < 0.5 & all_para_ck < 0.5) %>%
  filter(n <5) %>% # bottom 25%
  dplyr::count(dom1, dom2) %>%
  mutate(total = sum(n)) %>%
  mutate(pct_identity = n / total) %>%
  unite(col = dom1_dom2, c(dom1, dom2), sep = "_")
rin$dom1_dom2 <- factor(rin$dom1_dom2, levels = rin$dom1_dom2[order(rin$pct_identity)])
rin %>%
  ggplot(aes(dom1_dom2, pct_identity)) +
  geom_bar(stat = "identity") +
  # geom_line(aes(x = dom, y = 1-pct_identity, group = 1)) +
  # geom_point(shape = 21, colour = "white", fill = "black", size = 5, stroke = 2) +
  # ylim(0,0.5)+
  coord_flip() +
  theme_minimal()

# ggsave(filename = "rewiring_paired.pdf",
#        plot = last_plot(), path = "output/F4S1/",
#        width = 4,
#        height = 6)

# how many NTc, NTr, ECL2?
temp <- rin %>% 
  filter(grepl("NTc|NTr|ECL2", dom1_dom2)) %>%
  mutate(newsum = sum(n)) %>%
  dplyr::select(newsum, total) %>%
  unique() %>%
  mutate(newpct = newsum / total)
paste0("N-termini and ECL2 account for ", round(temp$newpct*100, ), " percent of rewired contacts")
