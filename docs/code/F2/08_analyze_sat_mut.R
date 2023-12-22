# To identify mutation sites with known minimal impact on CXCL12 activation
# (ie negative control), we used Wescott, Supp Table 3-5 and identified CXCR4 
# mutations with minimal effects on signaling as those with ≥95% efficacy; only 
# residues with GPCRdb positions were selected; the list includes...
# 3.53, 3.54, 3.57, 5.64, 6.32, 8.47, 8.49 (Supp Table 3; Gpro interface)
# 5.36(x37), 5.37(x38), 5.40(x41), 5.44(x45) (Supp Table 4; dimerization interface)
# 3.51, 5.59, 5.62, 5.63, 5.64, 6.39, 6.43, 7.49, 7.52 (Supp Table 5; signaling motifs)
# Of these residues, the TM5 ones are at the interface region but do not
# make structural contacts (contacts are 5x33, 5x36, 5x39, 5x40)
#
# As "positive control", from Stephens 2020 Sci Sig paper, 6x58 mutations have 
# most profound effects on CXCL12 binding and are at the interface

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import
data <- read_csv("data/mutagenesis_heredia/processed/cxcr4_clean_gpcrdb_means_log.csv")
data2 <- data
  
# shapiro wilk to test for normality
test <- data %>% filter(gn %in% c("5x38")) %>% 
  filter(sele == "cxcl12") %>%
  dplyr::select(resno, gn, value) %>% unique()
shapiro.test(test$value) # p-value = 0.001461 => p < 0.05 - NOT NORMALLY DISTRIBUTED

test <- data %>% filter(gn %in% c("6x58")) %>% 
  filter(sele == "cxcl12") %>%
  dplyr::select(resno, gn, value) %>% unique()
shapiro.test(test$value) # p-value = 0.004244 => p < 0.05 - NOT NORMALLY DISTRIBUTED

test <- data %>% filter(gn %in% c("NTr.Cm1")) %>% 
  filter(sele == "cxcl12") %>%
  dplyr::select(resno, gn, value) %>% unique()
shapiro.test(test$value) # p-value = 0.2034 

test <- data %>% filter(gn %in% c("1x22")) %>% 
  filter(sele == "cxcl12") %>%
  dplyr::select(resno, gn, value) %>% unique()
shapiro.test(test$value) # p-value = 0.6424 

test <- data %>% filter(gn %in% c("7x24")) %>% 
  filter(sele == "cxcl12") %>%
  dplyr::select(resno, gn, value) %>% unique()
shapiro.test(test$value) # p-value = 0.01934 => p < 0.05 - NOT NORMALLY DISTRIBUTED

# not normally distributed - use non-parametric test for pairwise comparisons
# https://www.datanovia.com/en/lessons/kruskal-wallis-test-in-r/
res.kruskal <- test %>% kruskal_test(value ~ gn)
res.kruskal
test %>% kruskal_effsize(value ~ gn)
#  Dunn's test instead of Wilcoxson for identifying pairwise p-values
# "Compared to the Wilcoxon’s test, the Dunn’s test takes into account the 
# rankings used by the Kruskal-Wallis test. It also does ties adjustments."
# Pairwise comparisons
pwc <- test %>% 
  dunn_test(value ~ gn, p.adjust.method = "bonferroni") 
pwc

# plot - BOXPLOT
data$gn <- factor(data$gn, levels = c("5x38", "6x58","NTr.Cm1", "1x22", "7x24"))
data %>% filter(gn %in% c("5x38", "6x58",  "NTr.Cm1", "1x22", "7x24")) %>%
  dplyr::select(resno, gn, sele, sub_mean) %>%
  filter(sele == "cxcl12") %>%
  unique() %>%
  ggplot(aes(gn, sub_mean*-1)) +
  geom_boxplot() +
  theme_minimal()

# ggsave(filename = "heredia_boxplot.pdf",
#        plot = last_plot(), path = "output/F2/",
#        width = 3,
#        height = 4.5)

