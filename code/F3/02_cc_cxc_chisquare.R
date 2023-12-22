# Takes CC and CXC chemokines/GPCRs with interaction strength >= 2 and 
# calculates statistical preference for intra-subfamily (e.g. CC-with-CC)
# interactions over inter-subfamily interactions

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# import, select those with evidence
data <- read_csv("data/network/positive_interactions.csv") %>%
  filter(grepl("CC", chemokine) | grepl("CXC", chemokine)) %>%
  filter(grepl("CC", gpcr) | grepl("CXC", gpcr))


# add class labels
test <- data %>% mutate(ck_class = case_when(
  grepl("CC", data$chemokine) ~ "cc_ck",
  grepl("CXC", data$chemokine) ~ "cxc_ck"
)) %>% mutate(ckr_class = case_when(
  grepl("CC", data$gpcr) ~ "cc_ckr",
  grepl("CXC", data$gpcr) ~ "cxc_ckr"
))

# count
test <- test %>% dplyr::select(ck_class, ckr_class) %>% dplyr::count(ck_class, ckr_class)
test <- test %>% pivot_wider(names_from = ck_class, values_from = n)
temp <- test$ckr_class
test <- test[,-1]
rownames(test) <- temp
rm(temp)

# stats - chi squared 
# http://www.sthda.com/english/wiki/chi-square-test-of-independence-in-r
test <- chisq.test(test)
test
test$observed
test$expected
test$residuals
library(corrplot)
corrplot(test$residuals, is.cor = FALSE)
contrib <- 100*test$residuals^2/test$statistic
contrib
corrplot(contrib, is.cor = FALSE)
test$p.value
paste0("The Chi-squared p-value is ", test$p.value)


