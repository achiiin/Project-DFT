# This is designed to find Machine Learning algorithm for predicting new compounds
# would be outliers or not.
# Class: 1 > outlier  0 > fitted 
# Input: 1. CountBB_0427.csv   (descirbe how SMILES is composed of building blocks) 
#        2. DFT_HOMO_BP86s_1percent_1outlier0fitted.csv (indicates which SMILES are outliers)

select_flavors <- "F18"
df_bb <- read.csv('BB_test_2.csv')
df_outlier <- read.csv('DFT_HOMO_BP86s_1percent_1outlier0fitted.csv')


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

inTrain = createDataPartition(df_bb$Class, p = 3/5)[[1]]
training = df_bb[ inTrain,-c(1)]
testing_1 = df_bb[-inTrain,-c(1)]
inTest = createDataPartition(testing_1$Class, p = 1/2)[[1]]
testing = testing_1[inTest,]
validation = testing_1[-inTest,]

####Random Forest####
fitControl <- trainControl(method = "none")
tgrid <- expand.grid(mtry=c(6)) 
Fit_rf <- train(Class~., trControl = fitControl, tuneGrid=tgrid,data=training,
                method = 'rf')

print(Fit_rf);print(Fit_rf$finalModel)

####glm####
Fit_glm <- train(Class~.,data=training,method = 'glm')

print(Fit_glm);print(Fit_glm$finalModel)






