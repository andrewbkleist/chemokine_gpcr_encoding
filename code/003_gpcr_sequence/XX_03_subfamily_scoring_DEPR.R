# Imports sequence matrix, partitions sequences into training and test sets,
# trains position-specific logistic regression models based on training data,
# and evaluates accuracy of trained model at accurately predicting subfamily
# (i.e. CC versus CXC) of the blinded sequence. Also uses the same position-
# specific, trained logistic regression models to assign prediction probability 
# scores (these are logistic regression output parameters used to make the 
# prediction) to all positions for human class A paralogs and for a subset
# of viral chemokine GPCR sequences. Note that the model can only evaluate test
# set sequences which contain residues for which it has been trained. To avoid 
# generating test sets that contained residues at a particular position 
# that were not present at the same position in the training data, random seeds 
# that partitioned the training/test sets accordingly were identified 
# empirically.
# 
# Note that the code below is a copy of code was run in a local directory on 
# 11.01.2020, and the output (test set accuracy and SD of logistic regression 
# models for each alignment position was copied over to 
# /data/sequence/gpcr/processed and not generated from this directory.
# This was done because the code utilizes functions from the package 
# InformationValue, which is no longer supported by R. The code below includes
# the original paths and also embeds functions within the code instead of 
# placing them within the centralized /code/000_functions.R script.

################################################################################

# ##### FUNCTIONS ################################################################
# 
# # FUNCTION 1 -----------------------------------------------------------------
# BuildTrainTest <- function(FILE, SEED, TRAIN.PCT){
#   
#   data <- read.csv(FILE, colClasses = "factor") # all colums *MUST* be factors (20200922 troubleshoot)
#   data <- data[, sapply(data, nlevels) > 1] # avoids errors from single factor columns
#   data <- data %>% select(-protein, -seq)
#   
#   # define classes
#   A <- data %>% filter(class == "cc") 
#   B <- data %>% filter(class == "cxc")
#   
#   # random selection of training rows
#   set.seed(SEED)  # for repeatability of samples
#   rows.A <- sample(1:nrow(A), TRAIN.PCT*nrow(B))  # deliberate B for even sampling
#   rows.B <- sample(1:nrow(B), TRAIN.PCT*nrow(B)) 
#   train.A <- A[rows.A, ]  
#   train.B <- B[rows.B, ]
#   train <- rbind(train.A, train.B)  # row bind the 1's and 0's 
#   
#   # define test rows
#   test.A <- A[-rows.A, ]
#   test.B <- B[-rows.B, ] 
#   test.A <- sample_n(test.A, nrow(test.B)) # match row numbers
#   
#   test <- rbind(test.A, test.B)  # row bind the 1's and 0's 
#   train$cat <- c("train")
#   test$cat <- c("test")
#   all <- rbind(train, test)
#   return(all)
#   rm(A, B, rows.A, rows.B, train.A, train.B, all, test.A, test.B, all)
# }
# 
# 
# # FUNCTION 2 -----------------------------------------------------------------
# LogResAcc <- function(TRAIN, TEST){
#   results <- data.frame()
#   
#   for(i in names(TRAIN)[2:ncol(TRAIN)]){
#     
#     # make model
#     fmla <- as.formula(paste0("class ~ ", i)) 
#     logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
#     
#     # **KEY LINES** OVERIDE FACTORS TO INCLUDE NEW FACTORS THAT ARE 
#     # ONLY IN TEST SET BUT NOT TRAINING  - 20200908
#     cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(as.factor(logitMod$xlevels$", i, "), as.factor(test$", i, "))"))
#     eval(parse(text = cmd))
#     
#     # evaluate model on test / convert to table
#     predicted <- plogis(predict(logitMod, test))  # predicted scores
#     con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics
#     acc <- (con.mat[1,1] + con.mat[2,2]) / 
#       (con.mat[1,1] + con.mat[2,2] + con.mat[1,2] + con.mat[2,1])
#     acc <- as.data.frame(acc)
#     acc$motif <- paste(i)
#     results <- rbind(results, acc)
#     rm(logitMod, predicted, con.mat, acc)
#     
#   }
#   results <- results %>% arrange(desc(acc))
#   return(results)
#   rm(results)
# }
# 
# # FUNCTION 3 -----------------------------------------------------------------
# LogResScore <- function(TRAIN){
#   results <- data.frame()
#   
#   test.master <- read.csv("03_ckr_seq/data/processed/ALL_classa_df.csv") 
#   #test.master <- subset(test.master, grepl("human", test.master$seq))
#   map.class <- test.master %>% select(protein, class)
#   test.master$class <- as.character(test.master$class)
#   protein.names <- test.master$protein
#   protein.names <- as.character(protein.names)
#   #test.master$class <- as.character(test.master$class)
#   
#   test.master <- test.master %>% mutate(class = case_when(
#     class == "non" ~ "cc",
#     class != "non" ~ class
#   ))
#   
#   for(j in protein.names){
#     
#     test <- test.master %>% filter(protein == j)
#     test <- test %>% select(-protein, -seq)
#     
#     for(i in names(TRAIN)[2:ncol(TRAIN)]){
#       
#       fmla <- as.formula(paste0("class ~ ", i)) 
#       logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
#       
#       # **KEY LINES** OVERIDE FACTORS TO INCLUDE NEW FACTORS ONLY IN TEST SET BUT NOT TRAINING  - 20200908
#       cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(as.factor(logitMod$xlevels$", i, "), as.factor(test$", i, "))"))
#       eval(parse(text = cmd))
#       
#       predicted <- plogis(predict(logitMod, test))  # predicted scores
#       con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics
#       
#       df <- as.data.frame(predicted)
#       df$position <- paste(i)
#       df$protein <- paste(j)
#       results <- rbind(results, df)
#       rm(logitMod, predicted, con.mat, df)
#     }
#     rm(test)
#   }
#   colnames(results)[1] <- c("cc_cxc")
#   results$class <- map.class$class[match(unlist(results$protein), map.class$protein)]
#   return(results)
#   rm(results)
# }
# 
# # FUNCTION 4 -----------------------------------------------------------------
# LogResScoreVirus <- function(TRAIN){
#   results <- data.frame()
#   
#   test.master <- read.csv("03_ckr_seq/data/processed/ALL_virus_df.csv") 
#   #test.master <- subset(test.master, grepl("human", test.master$seq))
#   map.class <- test.master %>% select(protein, class)
#   test.master$class <- as.character(test.master$class)
#   protein.names <- test.master$protein
#   protein.names <- as.character(protein.names)
#   #test.master$class <- as.character(test.master$class)
#   
#   test.master <- test.master %>% mutate(class = case_when(
#     class == "non" ~ "cc",
#     class != "non" ~ class
#   ))
#   
#   for(j in protein.names){
#     
#     test <- test.master %>% filter(protein == j)
#     test <- test %>% select(-protein, -seq)
#     
#     for(i in names(TRAIN)[2:ncol(TRAIN)]){
#       
#       fmla <- as.formula(paste0("class ~ ", i)) 
#       logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
#       
#       # **KEY LINES** OVERIDE FACTORS TO INCLUDE NEW FACTORS ONLY IN TEST SET BUT NOT TRAINING  - 20200908
#       cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(as.factor(logitMod$xlevels$", i, "), as.factor(test$", i, "))"))
#       eval(parse(text = cmd))
#       
#       predicted <- plogis(predict(logitMod, test))  # predicted scores
#       con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics
#       
#       df <- as.data.frame(predicted)
#       df$position <- paste(i)
#       df$protein <- paste(j)
#       results <- rbind(results, df)
#       rm(logitMod, predicted, con.mat, df)
#     }
#     rm(test)
#   }
#   colnames(results)[1] <- c("cc_cxc")
#   results$class <- map.class$class[match(unlist(results$protein), map.class$protein)]
#   return(results)
#   rm(results)
# }
# 
# 
# ##### 1: SCORE CKR SEQUENCES ###################################################
# 
# # RUN 1
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 100, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# test <- data %>% filter(cat == "test") %>% select(-cat)
# accLR.ckr1 <- LogResAcc(train, test)
# rm(data, train, test) 
# 
# # RUN 2
# # import, split to training, test sets
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 200, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# test <- data %>% filter(cat == "test") %>% select(-cat)
# accLR.ckr2 <- LogResAcc(train, test)
# rm(data, train, test)
# 
# # RUN 3
# # import, split to training, test sets
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 300, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# test <- data %>% filter(cat == "test") %>% select(-cat)
# accLR.ckr3 <- LogResAcc(train, test)
# rm(data, train, test)
# 
# # COMBINE, STD DEV
# colnames(accLR.ckr1)[1] <- c("acc1")
# colnames(accLR.ckr2)[1] <- c("acc2")
# colnames(accLR.ckr3)[1] <- c("acc3")
# 
# master <- left_join(accLR.ckr1, accLR.ckr2)
# master <- left_join(master, accLR.ckr3)
# master <- master %>% dplyr::select(motif, acc1, acc2, acc3)
# master <- master %>% gather(acc, value, 2:4)
# master <- master %>% group_by(motif) %>% summarise(mean = mean(value), sd = sd(value))
# master <- master %>% filter(!is.na(mean))
# 
# # write output
# write_csv(master, "03_ckr_seq/output/CKR_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# # **NOTE** ONLY COLS FOR WHICH >1 VARIABLE ARE EXIST ARE INCLUDED,
# # AS ARE COLS FOR WHICH ALL 3 REPLICATES GAVE VALUES
# rm(accLR.ckr1, accLR.ckr2, accLR.ckr3, master, BuildTrainTest, LogResAcc)
# 
# 
# ##### 2: PER RECEPTOR PER POSITION SCORES (ALL CLASS A INCLUDING CKRS) #########
# 
# # RUN 1
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 80, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.1 <- LogResScore(train)
# scoreLR.1$run <- c("run1")
# 
# # RUN 2
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 90, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.2 <- LogResScore(train)
# scoreLR.2$run <- c("run2")
# 
# # RUN 3
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 100, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.3 <- LogResScore(train)
# scoreLR.3$run <- c("run3")
# 
# rm(train)
# 
# # COMBINE
# scr <- bind_rows(scoreLR.1, scoreLR.2, scoreLR.3)
# scr <- scr %>% group_by(position, protein) %>% mutate(mean = mean(cc_cxc), sd = sd(cc_cxc)) %>% ungroup()
# 
# write_csv(scr, "03_ckr_seq/output/CKR_CLASSA_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# # rm(scoreLR.1, scoreLR.2, scoreLR.3)
# # rm(train, data)
# 
# ##### 3: PER RECEPTOR PER POSITION SCORES (VIRUS) ######################################
# 
# # RUN 1
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 100, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.1 <- LogResScoreVirus(train)
# scoreLR.1$run <- c("run1")
# 
# # RUN 2
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 200, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.2 <- LogResScoreVirus(train)
# scoreLR.2$run <- c("run2")
# 
# # RUN 3
# # import, select training, test, run
# data <- BuildTrainTest("03_ckr_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 300, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.3 <- LogResScoreVirus(train)
# scoreLR.3$run <- c("run3")
# 
# rm(train)
# 
# # COMBINE
# scr <- bind_rows(scoreLR.1, scoreLR.2, scoreLR.3)
# scr <- scr %>% group_by(position, protein) %>% mutate(mean = mean(cc_cxc), sd = sd(cc_cxc)) %>% ungroup()
# scr <- scr %>% select(protein,class,position,run, cc_cxc,mean,sd)
# colnames(scr)[5] <- c("score")
# 
# write_csv(scr, "03_ckr_seq/output/CKR_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# # rm(scoreLR.1, scoreLR.2, scoreLR.3)
# # rm(train, data)
