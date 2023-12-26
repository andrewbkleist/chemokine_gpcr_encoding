# Uses   position-specific, trained logistic regression models to assign 
# prediction probability scores (these are logistic regression output parameters 
# used to make the prediction) to all positions for human class A GPCRs.
#
# The model can only evaluate test set sequences that which contain residues for 
# which it has been trained. To avoid generating test sets that contained
# residues at a particular position that were not present at the same position 
# in the training data, test set rows (sequences) containing amino acids
# that were not present in training data were removed from the test set.
# Given random and equal partitioning of CC and CXC sequences to both
# sequence sets, test set pruning should not affect accuracy and should
# predominantly affect extreme N-termini.

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# RUN 1
# import, select training, test, run
data <- BuildTrainTest("data/sequence/gpcr/alignment_csv/ALL_cc_cxc_ortho_df.csv", 100, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat)
scoreLR.1 <- LogResScoreParalogClassA(train)
scoreLR.1$run <- c("run1")

# RUN 2
# import, select training, test, run
data <- BuildTrainTest("data/sequence/gpcr/alignment_csv/ALL_cc_cxc_ortho_df.csv", 200, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat)
scoreLR.2 <- LogResScoreParalogClassA(train)
scoreLR.2$run <- c("run2")

# RUN 3
# import, select training, test, run
data <- BuildTrainTest("data/sequence/gpcr/alignment_csv/ALL_cc_cxc_ortho_df.csv", 300, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat)
scoreLR.3 <- LogResScoreParalogClassA(train)
scoreLR.3$run <- c("run3")

rm(train)

# COMBINE
scr <- bind_rows(scoreLR.1, scoreLR.2, scoreLR.3)
scr <- scr %>% group_by(position, protein) %>% mutate(mean = mean(cc_cxc), sd = sd(cc_cxc)) %>% ungroup()

write_csv(scr, "03_ckr_seq/output/CKR_CLASSA_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# rm(scoreLR.1, scoreLR.2, scoreLR.3)
# rm(train, data)


