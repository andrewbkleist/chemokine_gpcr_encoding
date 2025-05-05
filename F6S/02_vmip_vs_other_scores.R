# packages, working directsory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

# import
data <- read_csv("02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
vmipii <- read_csv("02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
data <- rbind(data, vmipii)
t2 <- read_csv("02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv") %>%
  filter(mean >= 0.75)
t2 <- as.character(t2$motif)
rm(vmipii)

# order
data$position <- factor(data$position, levels = c("NTc.Cm4", "NTc.Cm1"))

# plot
data %>%
  filter(position %in% c("NTc.Cm1", "NTc.Cm4")) %>%
  filter(protein %in% c("ccl8", "vmip2xhhv8p")) %>%
  select(protein, position, mean, sd) %>%
  
  ggplot(aes(position, mean, fill = protein)) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd, width=.2)) +
  geom_point(shape = 21,   size = 4, stroke = 0.5) +
  theme_minimal() +
  # ylim(0,1) +
  # scale_fill_manual(values=c("steelblue4", "red4")) +
  coord_flip() +
  theme(axis.text.x=element_text(angle=90,vjust=.2, hjust=0))

ggsave(filename = "mutation_score_change.pdf", 
       plot = last_plot(), path = "F6S/output/",
       width = 4,
       height = 3)

  
  