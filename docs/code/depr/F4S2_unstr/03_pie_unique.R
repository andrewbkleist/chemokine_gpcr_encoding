source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/structure/processed/RIN_pairwise.csv") %>%
  filter(type %in% c("a_not_b","b_not_a")) 

data <- data %>%
  mutate(str_unstr_ck = case_when(
    grepl("NTc",source_gnccn) ~ "unstr",
    !grepl("NTc",source_gnccn) ~ "str")) %>%
  mutate(str_unstr_ckr = case_when(
    grepl("NTr",target_gnccn) ~ "unstr",
    !grepl("NTr",target_gnccn) ~ "str")) %>%
    dplyr::count(str_unstr_ck, str_unstr_ckr)   

data <- data %>%
  unite(col = type, c(str_unstr_ck, str_unstr_ckr), sep = "_")

data <- data %>% mutate(pct = n / sum(n))

hsize <- 1
data %>%
  ggplot(aes(x = 1, y= pct, fill = type)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()
  

# ---
data <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>%
  filter(!(file %in% c("6meo"))) %>%
  filter(no_pdb == 1 ) %>%
  select(source_gnccn, target_gnccn, dom1, dom2) %>%
  unique()

data <- data %>%
  mutate(str_unstr_ck = case_when(
    dom1 == "NTc" ~ "unstr",
    dom1 != "NTc" ~ "str",
  )) %>% mutate(str_unstr_ckr = case_when(
    dom2 == "NTr" ~ "unstr",
    dom2 != "NTr" ~ "str",
  )) %>% 
  dplyr::count(str_unstr_ck, str_unstr_ckr) 

data <- data %>%
  unite(col = type, c(str_unstr_ck, str_unstr_ckr), sep = "_")

data <- data %>% mutate(pct = n / sum(n))

data %>%
  ggplot(aes(x = "", y= pct)) +
  geom_bar(width = 0.5,size = 0.5, stat="identity", color = "white") +
  coord_polar("y") +
  theme_classic() +
  theme(axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

# ggsave(filename = "str_unstr_pie.pdf",
#        plot = last_plot(), path = "output/F4S2/",
#        width = 3,
#        height = 3)
  