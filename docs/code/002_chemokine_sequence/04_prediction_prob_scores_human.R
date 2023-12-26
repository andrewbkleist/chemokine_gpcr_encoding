# Uses   position-specific, trained logistic regression models to assign 
# prediction probability scores (these are logistic regression output parameters 
# used to make the prediction) to all positions for human chemokine paralogs.
#
# The model can only evaluate test set sequences that which contain residues for 
# which it has been trained. To avoid generating test sets that contained
# residues at a particular position that were not present at the same position 
# in the training data, test set rows (sequences) containing amino acids
# that were not present in training data were removed from the test set.
# Given random and equal partitioning of CC and CXC sequences to both
# sequence sets, test set pruning should not affect accuracy and should
# predominantly affect extreme N-termini.
#
# Random seeds adjusted in some instances to avoid the following error:
# Error in `contrasts<-`(`*tmp*`, value = contr.funs[1 + isOF[nn]]) : 
# contrasts can be applied only to factors with 2 or more levels
# This is a consequence of training/test set partitions which fail to deliver
# at least 2 distinct amino acids (including "-") and thus the model
# is not trainable by glm. Partitions adjusted empirically and blinded to
# results. This error mostly affects very highly gapped positions at extreme
# termini.
# https://stackoverflow.com/questions/18171246/error-in-contrasts-when-defining-a-linear-model-in-r

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

# RUN 1
# import, select training, test, run
data <- BuildTrainTest("data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv", 60, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat) # %>% dplyr::select(class, NTc.Cm1:CX.1)
scoreLR.1 <- LogResScoreParalogCK(train)
scoreLR.1$run <- c("run1")

# RUN 2
# import, select training, test, run
data <- BuildTrainTest("data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv", 70, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat)
scoreLR.2 <- LogResScoreParalogCK(train)
scoreLR.2$run <- c("run2")

# RUN 3
# import, select training, test, run
data <- BuildTrainTest("data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv", 81, 0.8)
train <- data %>% filter(cat == "train") %>% select(-cat)
scoreLR.3 <- LogResScoreParalogCK(train)
scoreLR.3$run <- c("run3")

rm(train)

# COMBINE
scr <- bind_rows(scoreLR.1, scoreLR.2, scoreLR.3)
scr <- scr %>% group_by(position, protein) %>% mutate(mean = mean(cc_cxc), sd = sd(cc_cxc)) %>% ungroup()
scr <- scr %>% dplyr::select(protein,class,position,run, cc_cxc,mean,sd)
colnames(scr)[5] <- c("score")

# write_csv(scr, "data/sequence/chemokine/processed/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv") # 20231223
rm(scoreLR.1, scoreLR.2, scoreLR.3)
rm(scr, data)

