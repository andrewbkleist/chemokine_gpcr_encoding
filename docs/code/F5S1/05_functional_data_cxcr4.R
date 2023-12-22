source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################
  
  data <- read_csv("data/functional/from_other_papers/ziarek_2013/ziarek_2013.csv") %>% 
    dplyr::select(variant, max_rat, ec50_rat)
  
  order <- c("WT", "Y7A", "Y12A", "Y21A")
  data$variant <- factor(data$variant, levels = rev(order))
  
  data %>%
    ggplot(aes(variant, ec50_rat)) +
    geom_bar(stat="identity") +
    theme_minimal() +
    coord_flip()
  
  # ggsave(filename = "ziarek.pdf",
  #        plot = last_plot(), path = "output/F5S1/",
  #        width = 3,
  #        height = 4)
  