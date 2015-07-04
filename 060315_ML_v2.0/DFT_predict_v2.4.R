# This is designed to find Machine Learning algorithm for predicting new compounds
# would be outliers or not.
# Class: 1 > outlier  0 > fitted 
# Input: 1. CountBB_BB_0508.csv   (descirbe how SMILES is composed of building blocks) 
#        2. --- (indicates which SMILES are outliers)

# v2.2.: use DFT_HOMO_BP86s_multi_1outlier0fitted_F18_V3.2.csv as input file.
#        only consider class a,d,e,f, regroup a as "0", and d,e,f as "1"
# v2.3. : add the measurements such as the area under the ROC curve, the sensitivity,
#         and spcificity
# v2.4. : Optimizing probability thresholds for class imbalances
#         (http://topepo.github.io/caret/custom_models.html)

select_flavors <- "F18"
# df_bb <- read.csv('BB_test_3.csv')
# df_outlier <- read.csv('BB_test_4.csv')
df_bb <- read.csv('CountBB_BB_0508.csv')
df_outlier <- read.csv('DFT_HOMO_BP86s_multi_1outlier0fitted_F18_V3.2.csv')


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
df_bb_2$Class[df_bb_2$Class %in% "a"] <- "0"
df_bb_2$Class[df_bb_2$Class %in% c("d","e","f")] <- "1"
df_bb_2$Class <- as.factor(df_bb_2$Class)
print(summary(df_bb_2$Class))

inTrain_80 = createDataPartition(df_bb_2$Class, p = 4/5)[[1]]
training_80 = df_bb_2[ inTrain_80,-c(1)]
testing_80 = df_bb_2[-inTrain_80,-c(1)]
print(names(training_80))
print(summary(training_80$Class))

##############################################






