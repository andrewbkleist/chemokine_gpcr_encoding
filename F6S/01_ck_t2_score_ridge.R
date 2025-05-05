# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)
library(ggridges)

##### 2: SCORE COMPARISONS ALL POSITIONS - RIDGE PLOT ##########################
  
# import
data <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
vmipii <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
data <- rbind(data, vmipii)
t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  filter(mean >= 0.75)
t2 <- as.character(t2$motif)
rm(vmipii)

# order
order <- c("ccl2", "ccl3", "ccl5", "ccl15", "ccl20", "vmip2xhhv8p",  "cxcl12", "cxcl8")
data$protein <- factor(data$protein, levels = rev(order))

# plot
data %>%
  filter(position %in% t2) %>%
  filter(protein %in% c("ccl2", "ccl3", "ccl5", "ccl15", "ccl20","vmip2xhhv8p",  "cxcl12", "cxcl8")) %>%
  ggplot(aes(x = mean, y = protein)) +
  geom_density_ridges(fill = "steelblue3") +
  scale_x_continuous(breaks=c(0, 0.5, 1)) +
  theme_minimal() 

ggsave(filename = "ridge_ck.pdf", 
       plot = last_plot(), path = "F6S/output/",
       width = 4,
       height = 5)

