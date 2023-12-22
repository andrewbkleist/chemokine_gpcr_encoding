source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import
data <- read_csv("data/network/Supplementary_Table_1.csv") %>%
  dplyr::select(chemokine, gpcr, author_date_journal, PMID, binding_kd_ki_nm, 
                binding_ec50_ic50_nm, signaling_ec50_ic50_nm, chemotaxis_ec50_ic50_nm, 
                chemotaxis_max_nm, ligand_type, interaction_strength, evidence_grade)


# pairings considered
ck.gpcr.pairs <- data %>% dplyr::select(chemokine, gpcr) %>% unique()
pct <- data_frame(name = c("pairs_considered", "pairs_not_considered"),
                  value = c( (nrow(ck.gpcr.pairs)/(46*23)), 1-(nrow(ck.gpcr.pairs)/(46*23)) ))

# make pie chart
hsize <- 1
p <- pct %>% 
  ggplot(aes(x = 1, y = value, fill = name)) +
  geom_col() +
  coord_polar(theta = "y") +
  xlim(c(0.2, hsize + 0.5)) +
  theme_minimal()
p

# ggsave(filename = "network_pie.pdf", 
#        plot = last_plot(), path = "output/F1S/",
#        width = 5,
#        height = 4)
