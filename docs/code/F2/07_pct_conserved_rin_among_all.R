source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import
data <- read_csv("data/integrated/RIN_CONS_CLASS.csv") %>% 
  filter(file != "6meo")
data <- data %>% unite(col = "gnccn", source_gnccn, target_gnccn, sep = "_", remove = FALSE)

# define conserved
temp <- data %>% 
  filter((all_para_ck >= 0.5) & (all_non_ackr_para_ckr >= 0.5))  %>%
  dplyr::select(gnccn) %>% 
  unique()
temp <- temp$gnccn
data <- data %>% 
  filter(file != "6meo") %>% 
  filter( !(is.na(all_non_ackr_para_ckr)) )  %>%
  filter(!(is.na(all_para_ck)) ) %>%

  mutate(conserved_rin = case_when(
  gnccn %in% temp ~ "yes",
  !(gnccn %in% temp) ~ "no",
))

rm(temp)

test <- data %>%
  dplyr::count(conserved_rin)

test$total <- nrow(data)
test <- test %>% mutate(pct = n / total)

# plot
hsize <- 1
test %>%
  ggplot(aes(x = 1, y = pct, fill = conserved_rin)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()

# ggsave(filename = "pct_rins_cons_pie.pdf", 
#        plot = last_plot(), path = "output/F2/",
#        width = 4,
#        height = 3)
