# packages, working directory
setwd("~/Desktop/Volkman_Lab/BABU_VOLKMAN/CHEMOKINE_PRJ/Data/chemokine_gpcr_ms/")
library(tidyverse)

################################################################################

data <- read_csv("08_ckr_motif/data/processed/20210710_motif_conversion_3mer_ckr_nterm.csv")
  
data <- subset(data, grepl(("DY|YD|EY|YE|D[A-Z]Y|Y[A-Z]D|E[A-Z]Y|Y[A-Z]E"), data$motif))

# select MOTIF "SLICE"
# data <- data %>% filter(motif %in% c("DYD", "EYD", "EYE", "DYE"))
  
# order.mot <- as.factor(unique(c(
#                                 "DY","YD", "EY", "YE",
#                                   "DxY","YxD", "ExY", "YxE",
#                                   "DxxY", "YxxD", "ExxY", "YxxE")))
# data$motif <- factor(data$motif, levels = rev(order.mot))

# order motifs based on appearance:
# want most unique to most conserved (global) then most to least conserved among orthologs
data <- data %>% add_count(motif)
data <- data %>% arrange(n, protein)
order.motif <- unique(data$motif)
levels(data$motif)

data$motif <- factor(data$motif, levels = rev(order.motif))

data$gn <- factor(data$gn, levels = rev(c("gnNTr.Cm28","gnNTr.Cm27","gnNTr.Cm26","gnNTr.Cm25","gnNTr.Cm24","gnNTr.Cm23","gnNTr.Cm22","gnNTr.Cm21","gnNTr.Cm20","gnNTr.Cm19","gnNTr.Cm18","gnNTr.Cm17","gnNTr.Cm16","gnNTr.Cm15","gnNTr.Cm14","gnNTr.Cm13","gnNTr.Cm12","gnNTr.Cm11","gnNTr.Cm10","gnNTr.Cm9","gnNTr.Cm8","gnNTr.Cm7","gnNTr.Cm6","gnNTr.Cm5","gnNTr.Cm4","gnNTr.Cm3","gnNTr.Cm2","gnNTr.Cm1")))

# order receptors
order.ckr <- as.factor(toupper(c("ccr1","ccr2", "ccr3", "ccr4",
                                 "ccr5","ccr6","ccr7","ccr8","ccr9",
                                 "ccr10","cxcr1","cxcr2","cxcr3",
                                 "cxcr4","cxcr5","cxcr6","cx3cr1", "xxcr1",
                                 "ackr1","ackr2","ackr3","ackr4",
                                 "ccrl2")))
levels(data$protein)
data$protein <- factor(data$protein, levels = (order.ckr))

data %>%
  ggplot(aes(gn, protein)) +
  geom_tile( fill = "mediumorchid4") +
  theme_minimal() +
  theme(axis.text.x=element_text(angle=-90,vjust=.2, hjust=0)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  coord_flip()

# ggsave(filename = "styr_motif_register.pdf", 
#        plot = last_plot(), path = "F5/output/",
#        width = 4,
#        height = 4)


