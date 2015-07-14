# This is designed to find Machine Learning algorithm for predicting new compounds
# would be outliers or not.
# Class: 1 > outlier  0 > fitted 
# Input: 1. CountBB_BB_0508.csv   (descirbe how SMILES is composed of building blocks) 
#        2. --- (indicates which SMILES are outliers)

# v2.2.: use DFT_HOMO_BP86s_multi_1outlier0fitted_F18_V3.2.csv as input file.
#        only consider class a,d,e,f, regroup a as "0", and d,e,f as "1"
# v2.3. : add the measurements such as the area under the ROC curve, the sensitivity,
#         and spcificity


select_flavors <- "F18"
df_bb <- read.csv('BB_test_3.csv')
df_outlier <- read.csv('BB_test_3_2.csv')
# df_bb <- read.csv('CountBB_BB_0508.csv')
# df_outlier <- read.csv('DFT_HOMO_BP86s_multi_1outlier0fitted_F18_V3.2.csv')


## check for missing packages and install them
if (!require("pacman")) install.packages("pacman")
pacman::p_load("caret", "rattle", "rpart.plot","gbm","e1071","nnet","pROC","MASS")
##
library(caret)
library(rattle)
library(rpart.plot)
library(gbm)
library(e1071)
library(nnet)
library(pROC)
library(MASS)


# list_1 <- df_bb$SMILES == df_outlier$SMILES
# 
# if(! FALSE %in% list_1){
#     df_bb <- data.frame(df_bb,Class = as.factor(df_outlier[,select_flavors]))
# }else{stop("Two files have different SMILES!")}

#####################################################
#### data preparation ####
##merging & process##
df_bb  <- merge(df_bb,df_outlier, by = "SMILES")
print(names(df_bb))

df_bb <- data.frame(df_bb,Class = as.factor(df_bb[,select_flavors]))
print(names(df_bb))

df_bb <- df_bb[ , -which(names(df_bb) %in% c("X",select_flavors))]
print(names(df_bb))
print(dim(df_bb))
###
set.seed(1990)

bad <- is.na(df_bb$Class)
df_bb <- df_bb[!bad,]
print(dim(df_bb))
print(summary(df_bb$Class))
inTrain = createDataPartition(df_bb$Class, p = 3/5)[[1]]
training = df_bb[ inTrain,-c(1)]
testing_1 = df_bb[-inTrain,-c(1)]
print(names(training))
print(names(testing_1))
set.seed(1990)
inTest = createDataPartition(testing_1$Class, p = 1/2)[[1]]
testing = testing_1[inTest,]
validation = testing_1[-inTest,]

#### create training(80%) & testing(20%) ####
# set.seed(1990)
# inTrain_80 = createDataPartition(df_bb$Class, p = 4/5)[[1]]
# training_80 = df_bb[ inTrain_80,-c(1)]
# testing_80 = df_bb[-inTrain_80,-c(1)]
# print(names(training_80))
# print(summary(training_80$Class))

#### regroup a as "0", and d,e,f as "1" ####
set.seed(1990)
df_bb_2 <- df_bb[df_bb$Class %in% c("a","d","e","f"),] 
print(summary(df_bb_2$Class))

df_bb_2$Class <- as.character(df_bb_2$Class)
df_bb_2$Class[df_bb_2$Class %in% "a"] <- "fit"
df_bb_2$Class[df_bb_2$Class %in% c("d","e","f")] <- "out"
df_bb_2$Class <- as.factor(df_bb_2$Class)
print(summary(df_bb_2$Class))

inTrain_80 = createDataPartition(df_bb_2$Class, p = 4/5)[[1]]
training_80 = df_bb_2[ inTrain_80,-c(1)]
testing_80 = df_bb_2[-inTrain_80,-c(1)]
print(names(training_80))
print(summary(training_80$Class))


#### Neural Networks ####
# set.seed(1990)
# Fit_nnet <- train(Class~., method="nnet", data = training_80,trace = FALSE)
# 
# print(Fit_nnet);print(Fit_nnet$finalModel)
# 
# pd_train_nnet <- predict(Fit_nnet,training_80)
# print(confusionMatrix(data = pd_train_nnet,reference = training_80$Class)$table)
# print(confusionMatrix(data = pd_train_nnet,reference = training_80$Class)$overall[1])
# 
# pd_test_nnet <- predict(Fit_nnet,testing_80)
# 
# print(confusionMatrix(data = pd_test_nnet,reference = testing_80$Class)$table)
# print(confusionMatrix(data = pd_test_nnet,reference = testing_80$Class)$overall[1])
####Random Forest####
set.seed(1990)
fitControl <- trainControl(method = "none",classProbs = TRUE,summaryFunction = twoClassSummary)
tgrid <- expand.grid(mtry=c(6)) 
Fit_rf <- train(Class~., trControl = fitControl, tuneGrid=tgrid,data=training_80,
                method = 'rf',ntree = 500, metric="ROC")

print(Fit_rf);print(Fit_rf$finalModel)

pd_test_rf <- predict(Fit_rf,testing_80)
print(confusionMatrix(data = pd_test_rf,reference = testing_80$Class)$table)
print(confusionMatrix(data = pd_test_rf,reference = testing_80$Class)$overall[1])
####glm####
# set.seed(1990)
# Fit_glm <- train(Class~.,data=training,method = 'glm')
# 
# print(Fit_glm);print(Fit_glm$finalModel)
# 
# pd_train_glm <- predict(Fit_glm,training)
# print(confusionMatrix(data = pd_train_glm,reference = training$Class)$table)
# print(confusionMatrix(data = pd_train_glm,reference = training$Class)$overall[1])
# 
# pd_test_glm <- predict(Fit_glm,testing)
# print(confusionMatrix(data = pd_test_glm,reference = testing$Class)$table)
# print(confusionMatrix(data = pd_test_glm,reference = testing$Class)$overall[1])
# ####SVM####
# set.seed(1990)
# # ctrl <- trainControl(method = "repeatedcv", number = 10, 
# #                      repeats = 5,savePred=T, classProb= T)
# # Fit_svm <- train(Class~.,data=training, method = "svmLinear", trControl = ctrl)
# Fit_svm <- svm(Class~.,data=training)
# print(Fit_svm);print(Fit_svm$finalModel)
# 
# pd_train_svm <- predict(Fit_svm,training)
# print(confusionMatrix(data = pd_train_svm,reference = training$Class)$table)
# print(confusionMatrix(data = pd_train_svm,reference = training$Class)$overall[1])
# 
# pd_test_svm <- predict(Fit_svm,testing)
# print(confusionMatrix(data = pd_test_svm,reference = testing$Class)$table)
# print(confusionMatrix(data = pd_test_svm,reference = testing$Class)$overall[1])


#### gbm ####
# set.seed(1990)
# Fit_gbm <- train(Class~., method="gbm",data = training)
# print(Fit_gbm);print(Fit_gbm$finalModel)
# 
# pd_train_gbm <- predict(Fit_gbm,training)
# print(confusionMatrix(data = pd_train_gbm,reference = training$Class)$table)
# print(confusionMatrix(data = pd_train_gbm,reference = training$Class)$overall[1])
# 
# pd_test_gbm <- predict(Fit_gbm,testing)
# print(confusionMatrix(data = pd_test_gbm,reference = testing$Class)$table)
# print(confusionMatrix(data = pd_test_gbm,reference = testing$Class)$overall[1])


#### lda ####
# set.seed(1990)
# 
# ctrl <- trainControl(method = "repeatedcv",
#                      repeats = 3,
#                      classProbs = TRUE,
#                      summaryFunction = twoClassSummary)
# 
# Fit_lda <- train(Class~.,data=training_80,method = 'lda',
#                  trControl = ctrl, metric="ROC")
# 
# print(Fit_lda);print(Fit_lda$finalModel)
# 
# pd_train_lda <- predict(Fit_lda,training_80)
# print(confusionMatrix(data = pd_train_lda,reference = training_80$Class)$table)
# print(confusionMatrix(data = pd_train_lda,reference = training_80$Class)$overall[1])
# 
# pd_test_lda <- predict(Fit_lda,testing_80)
# print(confusionMatrix(data = pd_test_lda,reference = testing_80$Class)$table)
# print(confusionMatrix(data = pd_test_lda,reference = testing_80$Class)$overall[1])




