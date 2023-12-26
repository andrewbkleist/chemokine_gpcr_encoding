# Imports sequence matrix, partitions sequences into training and test sets,
# trains position-specific logistic regression models based on training data,
# and evaluates accuracy of trained model at accurately predicting subfamily
# (i.e. CC versus CXC) of the blinded sequence. 
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
test <- data %>% filter(cat == "test") %>% dplyr::select(-cat) # %>% dplyr::select(class, NTc.Cm1:CX.1)
accLR.ckr1 <- LogResAcc(train, test)
rm(data, train, test)

# RUN 2
# import, split to training, test sets
data <- BuildTrainTest("data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv", 70, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat)
test <- data %>% filter(cat == "test") %>% dplyr::select(-cat)
accLR.ckr2 <- LogResAcc(train, test)
rm(data, train, test)
 
# RUN 3
# import, split to training, test sets
data <- BuildTrainTest("data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv", 81, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat)
test <- data %>% filter(cat == "test") %>% dplyr::select(-cat)
accLR.ckr3 <- LogResAcc(train, test)
rm(data, train, test)

# COMBINE, STD DEV
colnames(accLR.ckr1)[1] <- c("acc1")
colnames(accLR.ckr2)[1] <- c("acc2")
colnames(accLR.ckr3)[1] <- c("acc3")

master <- left_join(accLR.ckr1, accLR.ckr2)
master <- left_join(master, accLR.ckr3)
master <- master %>% dplyr::select(motif, acc1, acc2, acc3)
master <- master %>% pivot_longer(cols = 2:4, names_to = "repl")
master <- master %>% group_by(motif) %>% dplyr::summarise(mean = mean(value), sd = sd(value)) %>% ungroup()
master <- master %>% filter(!is.na(mean))

# write output
write_csv(master, "data/sequence/chemokine/processed/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# WRITTEN 2023.12.23
# **NOTE** ONLY COLS FOR WHICH >1 VARIABLE ARE EXIST ARE INCLUDED,
# AS ARE COLS FOR WHICH ALL 3 REPLICATES GAVE VALUES
rm(accLR.ckr1, accLR.ckr2, accLR.ckr3, master)

