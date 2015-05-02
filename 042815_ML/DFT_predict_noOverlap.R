# This is designed to find Machine Learning algorithm for predicting new compounds
# would be outliers or not.
# Class: 1 > outlier  0 > fitted 
# Input: 1. CountBB_0427.csv   (descirbe how SMILES is composed of building blocks) 
#        2. DFT_HOMO_BP86s_1percent_1outlier0fitted.csv (indicates which SMILES are outliers)

select_flavors <- "F18"
df_bb <- read.csv('BB_test_2.csv')
df_outlier <- read.csv('DFT_HOMO_BP86s_1percent_1outlier0fitted.csv')
# df_bb <- read.csv('countBB_0427_noOverlap.csv')
# df_outlier <- read.csv('DFT_HOMO_BP86s_1percent_1outlier0fitted_F18.csv')

## check for missing packages and install them
if (!require("pacman")) install.packages("pacman")
pacman::p_load("caret", "rattle", "rpart.plot")
##
library(caret)
library(rattle)
library(rpart.plot)



list_1 <- df_bb$SMILES == df_outlier$SMILES

if(! FALSE %in% list_1){
    df_bb <- data.frame(df_bb,Class = as.factor(df_outlier[,select_flavors]))
}else{stop("Two files have different SMILES!")}

#####################################################
#### data preparation ####
set.seed(1990)

bad <- is.na(df_bb$Class)
df_bb <- df_bb[!bad,]
df_bb_2 <- df_bb[,c(1,28)]
###
f_del0col = function(M) M[, colSums(abs(M)) != 0] ## for removing all columns with 0
df_bb <- f_del0col(df_bb[,-c(1,28)])
print(names(df_bb))
df_bb <- cbind(df_bb,df_bb_2)
print(names(df_bb))
###

inTrain = createDataPartition(df_bb$Class, p = 3/5)[[1]]
training = df_bb[ inTrain, !names(df_bb) %in% c("SMILES")]
print(dim(training))
print(names(training))
testing_1 = df_bb[-inTrain,!names(df_bb) %in% c("SMILES")]

set.seed(1990)
inTest = createDataPartition(testing_1$Class, p = 1/2)[[1]]
testing = testing_1[inTest,]
print(dim(testing))
print(names(testing))
validation = testing_1[-inTest,]
print(dim(validation))
print(names(validation))

####Random Forest####
# set.seed(1990)
# fitControl <- trainControl(method = "none")
# tgrid <- expand.grid(mtry=c(6)) 
# Fit_rf <- train(Class~., trControl = fitControl, tuneGrid=tgrid,data=training,
#                 method = 'rf',ntree = 500)
# 
# print(Fit_rf);print(Fit_rf$finalModel)

# ####glm####
set.seed(1990)
Fit_glm <- train(Class~.,data=training,method = 'glm')

print(Fit_glm);print(Fit_glm$finalModel)

pd_train_glm <- predict(Fit_glm,training)
print(confusionMatrix(data = pd_train_glm,reference = training$Class)$table)
print(confusionMatrix(data = pd_train_glm,reference = training$Class)$overall[1])

####SVM####
# set.seed(1990)
# ctrl <- trainControl(method = "repeatedcv", number = 10, 
#                      repeats = 5,savePred=T, classProb= T)
# Fit_svm <- train(Class~.,data=training, method = "svmLinear", trControl = ctrl)
# 
# print(Fit_svm);print(Fit_svm$finalModel)





