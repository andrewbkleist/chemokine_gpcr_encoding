# Imports sequence matrix, partitions sequences into training and test sets,
# trains position-specific logistic regression models based on training data,
# and evaluates accuracy of trained model at accurately predicting subfamily
# (i.e. CC versus CXC) of the blinded sequence. Also uses the same position-
# specific, trained logistic regression models to assign prediction probability 
# scores (these are logistic regression output parameters used to make the 
# prediction) to all positions for human chemokine paralogs and for a subset
# of viral chemokine sequences. Note that the model can only evaluate test
# set sequences which contain residues for which it has been trained. To avoid 
# generating test sets that contained residues at a particular position 
# that were not present at the same position in the training data, random seeds 
# that partitioned the training/test sets accordingly were identified 
# empirically.

source("code/000_libraries.R")
source("code/000_functions.R")

################################################################################

###### FUNCTIONS ################################################################

# FUNCTION 1 -----------------------------------------------------------------
BuildTrainTest <- function(FILE, SEED, TRAIN.PCT){

  data <- read.csv(FILE, colClasses = "factor") # all colums *MUST* be factors (20200922 troubleshoot)
  data <- data[, sapply(data, nlevels) > 1] # avoids errors from single factor columns by removing
  data <- data %>% dplyr::select(-protein, -seq)

  # define classes
  A <- data %>% filter(class == "cc")
  B <- data %>% filter(class == "cxc")

  # random selection of training rows
  set.seed(SEED)  # for repeatability of samples
  rows.A <- sample(1:nrow(A), TRAIN.PCT*nrow(B))  # deliberate B for even sampling
  rows.B <- sample(1:nrow(B), TRAIN.PCT*nrow(B))
  train.A <- A[rows.A, ]
  train.B <- B[rows.B, ]
  train <- rbind(train.A, train.B)  # row bind the 1's and 0's

  # define test rows
  test.A <- A[-rows.A, ]
  test.B <- B[-rows.B, ]
  test.A <- sample_n(test.A, nrow(test.B)) # match row numbers

  test <- rbind(test.A, test.B)  # row bind the 1's and 0's
  train$cat <- c("train")
  test$cat <- c("test")
  all <- rbind(train, test)
  return(all)
  rm(A, B, rows.A, rows.B, train.A, train.B, all, test.A, test.B, all)
}

 
# # FUNCTION 2 -----------------------------------------------------------------
LogResAcc <- function(TRAIN, TEST){
  results <- data.frame()

  for(i in names(TRAIN)[2:ncol(TRAIN)]){

    # (1) MAKE MODEL
    # see http://r-statistics.co/Logistic-Regression-With-R.html
      #--
      i <- "NTc.Cm69"
      TRAIN <- train
      TEST <- test
      #--
    fmla <- as.formula(paste0("class ~ ", i))
    logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))
    
    # (1.5) REMOVE ROWS FROM TEST CONTAINING FACTORS NOT IN TRAIN
    model.factors <- c(logitMod$xlevels$NTc.Cm69)
    TEST <- TEST %>% filter(NTc.Cm69 %in% c(model.factors))
    
    
    test2$NTc.Cm69 <- as.character(test2$NTc.Cm69) # convert to character
    levels(test2$NTc.Cm69) <- c(logitMod$xlevels$NTc.Cm69) # override
    test2$NTc.Cm69 <- as.factor(test2$NTc.Cm69) # comvert back to factor
    
    # (2) CHECK TRAINING & TEST FACTORS
    cmd <- noquote(paste0("factor.check <- setdiff(test$", i, ", logitMod$xlevels$", i, ")"))
    eval(parse(text = cmd))
    # factor.check <- setdiff(test$NTc.Cm69, logitMod$xlevels$NTc.Cm69)
    
    if(length(factor.check) > 0){
      acc <- NA
      acc <- as.data.frame(acc)
      acc$motif <- paste(i)
      results <- rbind(results, acc)
      rm(logitMod, acc, fmla)
      
    }
    if(length(factor.check) < 1){
      # evaluate model on test / convert to table
      predicted <- plogis(predict(logitMod, TEST)) # predicted scores
      # ...same as... predicted <- predict(logitMod, newdata = test, type = "response") 
      predicted <- factor(ifelse(predict(logitMod, TEST) > 0.5, "cxc", "cc"))
      #--
      # library(pROC)
      # plot.roc(test$class, predicted)
      #--
      
      # make confusion matrix, get accuracy
      con.mat <- caret::confusionMatrix(as.factor(predicted), as.factor(test$class))
      acc <- con.mat$overall[[1]] # accuracy; run single bracket to confirm: con.mat$overall[1]
      # sanity check:
      # temp <- cbind(data.frame(test_preds = predicted, test_actual = test$class))
      # temp <- temp %>% mutate(same = case_when( (test_preds == test_actual) ~ "yes")) 
      acc <- as.data.frame(acc)
      acc$motif <- paste(i)
      results <- rbind(results, acc)
      rm(logitMod, predicted, con.mat, acc, fmla)
    }
    
    # **KEY LINES** OVERIDE FACTORS TO INCLUDE NEW, TEST SET FACTORS
    # https://stackoverflow.com/questions/22315394/factor-has-new-levels-error-for-variable-im-not-using
    # cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(logitMod$xlevels$", i, ", levels(test$", i, "))"))
    # eval(parse(text = cmd))
      #--
      # logitMod$xlevels$NTc.Cm69 <- union(as.factor(logitMod$xlevels$NTc.Cm69), unique(as.character(test$NTc.Cm69)))
      # logitMod$xlevels$NTc.Cm69 <- union(logitMod$xlevels$NTc.Cm69, levels(test$NTc.Cm69))
      #--

  }
  results <- results %>% arrange(desc(acc))
  return(results)
  rm(results)
}

# # FUNCTION 3 -----------------------------------------------------------------
LogResScoreVirus <- function(TRAIN){
  results <- data.frame()

  test.master <- read.csv("02_ck_seq/data/processed/ALL_virus_df.csv")
  #test.master <- subset(test.master, grepl("human", test.master$seq))
  map.class <- test.master %>% select(protein, class)
  test.master$class <- as.character(test.master$class)
  protein.names <- test.master$protein
  protein.names <- as.character(protein.names)
  #test.master$class <- as.character(test.master$class)

  test.master <- test.master %>% mutate(class = case_when(
    class == "non" ~ "cc",
    class != "non" ~ class
  ))

  for(j in protein.names){

    test <- test.master %>% filter(protein == j)
    test <- test %>% select(-protein, -seq)

    for(i in names(TRAIN)[2:ncol(TRAIN)]){

      fmla <- as.formula(paste0("class ~ ", i))
      logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))

      # **KEY LINES** OVERIDE FACTORS TO INCLUDE NEW FACTORS ONLY IN TEST SET BUT NOT TRAINING  - 20200908
      cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(as.factor(logitMod$xlevels$", i, "), as.factor(test$", i, "))"))
      eval(parse(text = cmd))

      predicted <- plogis(predict(logitMod, test))  # predicted scores
      con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics

      df <- as.data.frame(predicted)
      df$position <- paste(i)
      df$protein <- paste(j)
      results <- rbind(results, df)
      rm(logitMod, predicted, con.mat, df)
    }
    rm(test)
  }
  colnames(results)[1] <- c("cc_cxc")
  results$class <- map.class$class[match(unlist(results$protein), map.class$protein)]
  return(results)
  rm(results)
}

# # FUNCTION 4 -----------------------------------------------------------------
LogResScoreParalog <- function(TRAIN){
  results <- data.frame()

  test.master <- read.csv("02_ck_seq/data/processed/ALL_para_df.csv")
  #test.master <- subset(test.master, grepl("human", test.master$seq))
  map.class <- test.master %>% select(protein, class)
  test.master$class <- as.character(test.master$class)
  protein.names <- test.master$protein
  protein.names <- as.character(protein.names)
  #test.master$class <- as.character(test.master$class)

  test.master <- test.master %>% mutate(class = case_when(
    class == "non" ~ "cc",
    class != "non" ~ class
  ))

  for(j in protein.names){

    test <- test.master %>% filter(protein == j)
    test <- test %>% select(-protein, -seq)

    for(i in names(TRAIN)[2:ncol(TRAIN)]){

      fmla <- as.formula(paste0("class ~ ", i))
      logitMod <- glm(fmla, data=TRAIN, family=binomial(link="logit"))

      # **KEY LINES** OVERIDE FACTORS TO INCLUDE NEW FACTORS ONLY IN TEST SET BUT NOT TRAINING  - 20200908
      cmd <- noquote(paste0("logitMod$xlevels$", i, " <- union(as.factor(logitMod$xlevels$", i, "), as.factor(test$", i, "))"))
      eval(parse(text = cmd))

      predicted <- plogis(predict(logitMod, test))  # predicted scores
      con.mat <- InformationValue::confusionMatrix(test$class, predicted) # see also ModelMetrics

      df <- as.data.frame(predicted)
      df$position <- paste(i)
      df$protein <- paste(j)
      results <- rbind(results, df)
      rm(logitMod, predicted, con.mat, df)
    }
    rm(test)
  }
  colnames(results)[1] <- c("cc_cxc")
  results$class <- map.class$class[match(unlist(results$protein), map.class$protein)]
  return(results)
  rm(results)
}


# ##### 1: SCORE CK SEQUENCES - CC VS. CXC ACCURACY ##############################
# 
# # RUN 1
# import, select training, test, run
data <- BuildTrainTest("data/sequence/chemokine/alignment_csv/ALL_cc_cxc_ortho_df.csv", 60, 0.8)
train <- data %>% filter(cat == "train") %>% dplyr::select(-cat) # %>% dplyr::select(class, NTc.Cm1:CX.1)
test <- data %>% filter(cat == "test") %>% dplyr::select(-cat) # %>% dplyr::select(class, NTc.Cm1:CX.1)
accLR.ckr1 <- LogResAcc(train, test)
rm(data, train, test)

# # RUN 2
# # import, split to training, test sets
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 79, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# test <- data %>% filter(cat == "test") %>% select(-cat)
# accLR.ckr2 <- LogResAcc(train, test)
# rm(data, train, test)
# 
# # RUN 3
# # import, split to training, test sets
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 101, 0.8)
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
# master <- master %>% group_by(motif) %>% summarise(mean = mean(value), sd = sd(value)) %>% ungroup()
# master <- master %>% filter(!is.na(mean))
# 
# # write output
# write_csv(master, "02_ck_seq/output/CK_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# # **NOTE** ONLY COLS FOR WHICH >1 VARIABLE ARE EXIST ARE INCLUDED,
# # AS ARE COLS FOR WHICH ALL 3 REPLICATES GAVE VALUES
# rm(accLR.ckr1, accLR.ckr2, accLR.ckr3, master)
# 
# 
# ##### 2: PER CHEMOKINE PER POSITION SCORES (CHEMOKINE HUMAN PARALOGS) ##########
# 
# # RUN 1
# # import, select training, test, run
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 60, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.1 <- LogResScoreParalog(train)
# scoreLR.1$run <- c("run1")
# 
# # RUN 2
# # import, select training, test, run
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 79, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.2 <- LogResScoreParalog(train)
# scoreLR.2$run <- c("run2")
# 
# # RUN 3
# # import, select training, test, run
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 101, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.3 <- LogResScoreParalog(train)
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
# write_csv(scr, "02_ck_seq/output/CK_PARALOG_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# rm(scoreLR.1, scoreLR.2, scoreLR.3)
# rm(scr, data)
# 
# ##### 3: PER CHEMOKINE PER POSITION SCORES (VIRUS) #############################
# 
# # RUN 1
# # import, select training, test, run
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 60, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.1 <- LogResScoreVirus(train)
# scoreLR.1$run <- c("run1")
# 
# # RUN 2
# # import, select training, test, run
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 79, 0.8)
# train <- data %>% filter(cat == "train") %>% select(-cat)
# scoreLR.2 <- LogResScoreVirus(train)
# scoreLR.2$run <- c("run2")
# 
# # RUN 3
# # import, select training, test, run
# data <- BuildTrainTest("02_ck_seq/data/processed/ALL_cc_cxc_ortho_df.csv", 101, 0.8)
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
# write_csv(scr, "02_ck_seq/output/CK_VIRUS_LOGISTIC_REGRESSION_ACCURACY_N3.csv")
# rm(scoreLR.1, scoreLR.2, scoreLR.3)
# rm(scr, data)
