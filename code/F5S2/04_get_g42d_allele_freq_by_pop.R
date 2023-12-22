source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

data <- read_csv("data/variant/gnomad/processed/CKR_GNOMAD_TABLE.csv") %>%
  filter(gene_symbol == "ACKR1" & Consequence == "p.Gly42Asp") %>%
  dplyr::select("Allele.Count.African", "Allele.Number.African",
                "Allele.Count.Ashkenazi.Jewish","Allele.Number.Ashkenazi.Jewish",
                "Allele.Count.East.Asian", "Allele.Number.East.Asian",
                "Allele.Count.European..Finnish.", "Allele.Number.European..Finnish.",
                "Allele.Count.Latino", "Allele.Number.Latino" , 
                "Allele.Count.South.Asian", "Allele.Number.South.Asian",
                "Allele.Count.Other", "Allele.Number.Other") 

data <- data %>%
  pivot_longer(cols = 1:ncol(data), names_to = c("a")) %>%
  mutate(no_vs_freq = case_when(
    grepl("Count", a) ~ "count",
    grepl("Number", a) ~ "no"
  ))
data <- data %>% 
  separate(col = a, into = c("a","b"), sep = "Allele.", remove = TRUE) %>%
  dplyr::select(-a) %>%
  separate(col = b, into = c("a","group"), sep = "\\.", remove = FALSE) %>%
  dplyr::select(-a, -b)
data <- data %>%
  pivot_wider(names_from = no_vs_freq, values_from = value) %>%
  mutate(pop_freq = count/no) %>%
  mutate(pop_freq_alt = 1-pop_freq) %>%
  pivot_longer(cols = c(pop_freq, pop_freq_alt), names_to  = "freq")

data$group <- factor(data$group, 
                     levels = rev(c("East", "South", "Latino", "European", 
                                "Ashkenazi","Other","African")))

data %>%
  ggplot(aes(fill = freq, x = value, y = group)) + # fill=condition, y=value, x=speci
  geom_bar(position="fill", stat="identity") +
  theme_minimal()
  
# ggsave(filename = "gnomad_ackr1_g42d_freq_pop.pdf",
#        plot = last_plot(), path = "output/F5S2/",
#        width = 7,
#        height = 3)
