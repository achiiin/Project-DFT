# This is designed to find Machine Learning algorithm for predicting new compounds
# would be outliers or not.
# Class: 1 > outlier  0 > fitted 
# Input: 1. CountBB_0427.csv   (descirbe how SMILES is composed of building blocks) 
#        2. DFT_HOMO_BP86s_1percent_1outlier0fitted.csv (indicates which SMILES are outliers)

select_flavors <- "F18"
# df_bb <- read.csv('BB_test_2.csv')
# df_outlier <- read.csv('DFT_HOMO_BP86s_1percent_1outlier0fitted.csv')
df_bb <- read.csv('CountBB_Existence_0510.csv')
df_outlier <- read.csv('DFT_HOMO_BP86s_1percent_1outlier0fitted_F18.csv')


## check for missing packages and install them
if (!require("pacman")) install.packages("pacman")
pacman::p_load("caret", "rattle", "rpart.plot","gbm","e1071")
##
library(caret)
library(rattle)
library(rpart.plot)
library(gbm)
library(e1071)



list_1 <- df_bb$SMILES == df_outlier$SMILES

if(! FALSE %in% list_1){
    df_bb <- data.frame(df_bb,Class = as.factor(df_outlier[,select_flavors]))
}else{stop("Two files have different SMILES!")}

#####################################################
#### data preparation ####
set.seed(1990)

bad <- is.na(df_bb$Class)
df_bb <- df_bb[!bad,]

inTrain = createDataPartition(df_bb$Class, p = 3/5)[[1]]
training = df_bb[ inTrain,-c(1)]
testing_1 = df_bb[-inTrain,-c(1)]

set.seed(1990)
inTest = createDataPartition(testing_1$Class, p = 1/2)[[1]]
testing = testing_1[inTest,]
validation = testing_1[-inTest,]

####Random Forest####
# set.seed(1990)
# fitControl <- trainControl(method = "none")
# tgrid <- expand.grid(mtry=c(6)) 
# Fit_rf <- train(Class~., trControl = fitControl, tuneGrid=tgrid,data=training,
#                 method = 'rf',ntree = 500)
# 
# print(Fit_rf);print(Fit_rf$finalModel)
# 
# pd_test_rf <- predict(Fit_rf,testing)
# print(confusionMatrix(data = pd_test_rf,reference = testing$Class)$table)
# print(confusionMatrix(data = pd_test_rf,reference = testing$Class)$overall[1])
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
####SVM####
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
# Fit_lda <- train(Class~.,data=training,method = 'lda')
# print(Fit_lda);print(Fit_lda$finalModel)
# 
# pd_train_lda <- predict(Fit_lda,training)
# print(confusionMatrix(data = pd_train_lda,reference = training$Class)$table)
# print(confusionMatrix(data = pd_train_lda,reference = training$Class)$overall[1])
# 
# pd_test_lda <- predict(Fit_lda,testing)
# print(confusionMatrix(data = pd_test_lda,reference = testing$Class)$table)
# print(confusionMatrix(data = pd_test_lda,reference = testing$Class)$overall[1])



