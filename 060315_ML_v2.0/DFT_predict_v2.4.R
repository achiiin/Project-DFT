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
# df_outlier <- read.csv('BB_test_3_2.csv')
df_bb <- read.csv('CountBB_BB_0508.csv')
df_outlier <- read.csv('DFT_HOMO_BP86s_multi_1outlier0fitted_F18_V3.2.csv')


## check for missing packages and install them
if (!require("pacman")) install.packages("pacman")
pacman::p_load("caret", "rattle", "rpart.plot","gbm","e1071","nnet","pROC",
               "MASS","plyr","reshape2")
##
library(caret)
library(rattle)
library(rpart.plot)
library(gbm)
library(e1071)
library(nnet)
library(pROC)
library(MASS)

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
# inTrain = createDataPartition(df_bb$Class, p = 3/5)[[1]]
# training = df_bb[ inTrain,-c(1)]
# testing_1 = df_bb[-inTrain,-c(1)]
# print(names(training))
# print(names(testing_1))
# set.seed(1990)
# inTest = createDataPartition(testing_1$Class, p = 1/2)[[1]]
# testing = testing_1[inTest,]
# validation = testing_1[-inTest,]

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
df_bb_2$Class[df_bb_2$Class %in% c("d","e","f")] <- "Out"
df_bb_2$Class <- as.factor(df_bb_2$Class)
print(summary(df_bb_2$Class))

# ### for 80% training & 20% testing ###
# inTrain_80 = createDataPartition(df_bb_2$Class, p = 4/5)[[1]]
# training_80 = df_bb_2[ inTrain_80,-c(1)]
# testing_80 = df_bb_2[-inTrain_80,-c(1)]
# print(names(training_80))
# print(dim(training_80))
# print(summary(training_80$Class))
# 
# ##############################################
# trainingSet <- training_80
# testingSet <- testing_80
# ###

### for 60% training & 20% testing & 20% validation###
inTrain_60 = createDataPartition(df_bb_2$Class, p = 3/5)[[1]]
training_60 = df_bb_2[ inTrain_60,-c(1)]
testing_40 = df_bb_2[-inTrain_60,-c(1)]
print(names(training_60))
print(dim(training_60))
print(summary(training_60$Class))
set.seed(1990)
inTest = createDataPartition(testing_40$Class, p = 1/2)[[1]]
testing_20 = testing_40[inTest,]
validation_20 = testing_40[-inTest,]
print(names(testing_20))
print(dim(testing_20))
print(summary(testing_20$Class))
##############################################
trainingSet <- training_60
testingSet <- testing_20
###

## Class frequencies
table(trainingSet$Class)

set.seed(949)

#################################
mod0 <- train(Class ~ ., data = trainingSet,
              method = "rf",
              metric = "ROC",
              tuneGrid = data.frame(mtry = 3),
              ntree = 500,
              trControl = trainControl(method = "repeatedcv",
                                       repeats = 5,
                                       classProbs = TRUE,
                                       summaryFunction = twoClassSummary))

# ###### for testing 
# mod0 <- train(Class ~ ., data = trainingSet,
#               method = "rf",
#               metric = "ROC",
#               tuneGrid = data.frame(mtry = 3),
#               ntree = 5,
#               trControl = trainControl(method = "repeatedcv",
#                                        repeats = 1,
#                                        classProbs = TRUE,
#                                        summaryFunction = twoClassSummary))
getTrainPerf(mod0)

## Get the ROC curve
roc0 <- roc(testingSet$Class,
            predict(mod0, testingSet, type = "prob")[,1],
            levels = rev(levels(testingSet$Class)))
roc0

## Now plot
png(filename = "ROC.png")
plot(roc0, print.thres = c(.5), type = "S",
     print.thres.pattern = "%.3f (Spec = %.2f, Sens = %.2f)",
     print.thres.cex = .8,
     legacy.axes = TRUE)
dev.off()

## Get the model code for the original random forest method:

thresh_code <- getModelInfo("rf", regex = FALSE)[[1]]
thresh_code$type <- c("Classification")
## Add the threshold as another tuning parameter
thresh_code$parameters <- data.frame(parameter = c("mtry", "threshold"),
                                     class = c("numeric", "numeric"),
                                     label = c("#Randomly Selected Predictors",
                                               "Probability Cutoff"))
## The default tuning grid code:
thresh_code$grid <- function(x, y, len = NULL) {
    p <- ncol(x)
    expand.grid(mtry = floor(sqrt(p)),
                threshold = seq(.5, .99, by = .01))
}

## Here we fit a single random forest model (with a fixed mtry)
## and loop over the threshold values to get predictions from the same
## randomForest model.
thresh_code$loop = function(grid) {
    library(plyr)
    loop <- ddply(grid, c("mtry"),
                  function(x) c(threshold = max(x$threshold)))
    submodels <- vector(mode = "list", length = nrow(loop))
    for(i in seq(along = loop$threshold)) {
        index <- which(grid$mtry == loop$mtry[i])
        cuts <- grid[index, "threshold"]
        submodels[[i]] <- data.frame(threshold = cuts[cuts != loop$threshold[i]])
    }
    list(loop = loop, submodels = submodels)
}

## Fit the model independent of the threshold parameter
thresh_code$fit = function(x, y, wts, param, lev, last, classProbs, ...) {
    if(length(levels(y)) != 2)
        stop("This works only for 2-class problems")
    randomForest(x, y, mtry = param$mtry, ...)
}

## Now get a probability prediction and use different thresholds to
## get the predicted class
thresh_code$predict = function(modelFit, newdata, submodels = NULL) {
    class1Prob <- predict(modelFit,
                          newdata,
                          type = "prob")[, modelFit$obsLevels[1]]
    ## Raise the threshold for class #1 and a higher level of
    ## evidence is needed to call it class 1 so it should 
    ## decrease sensitivity and increase specificity
    out <- ifelse(class1Prob >= modelFit$tuneValue$threshold,
                  modelFit$obsLevels[1],
                  modelFit$obsLevels[2])
    if(!is.null(submodels)) {
        tmp2 <- out
        out <- vector(mode = "list", length = length(submodels$threshold))
        out[[1]] <- tmp2
        for(i in seq(along = submodels$threshold)) {
            out[[i+1]] <- ifelse(class1Prob >= submodels$threshold[[i]],
                                 modelFit$obsLevels[1],
                                 modelFit$obsLevels[2])
        }
    }
    out
}

## The probabilities are always the same but we have to create
## mulitple versions of the probs to evaluate the data across
## thresholds
thresh_code$prob = function(modelFit, newdata, submodels = NULL) {
    out <- as.data.frame(predict(modelFit, newdata, type = "prob"))
    if(!is.null(submodels)) {
        probs <- out
        out <- vector(mode = "list", length = length(submodels$threshold)+1)
        out <- lapply(out, function(x) probs)
    }
    out
}


###################################

fourStats <- function (data, lev = levels(data$obs), model = NULL) {
    ## This code will get use the area under the ROC curve and the
    ## sensitivity and specificity values using the current candidate
    ## value of the probability threshold.
    out <- c(twoClassSummary(data, lev = levels(data$obs), model = NULL))
    
    ## The best possible model has sensitivity of 1 and specificity of 1. 
    ## How far are we from that value?
    coords <- matrix(c(1, 1, out["Spec"], out["Sens"]),
                     ncol = 2,
                     byrow = TRUE)
    colnames(coords) <- c("Spec", "Sens")
    rownames(coords) <- c("Best", "Current")
    c(out, Dist = dist(coords)[1])
}

set.seed(949)
mod1 <- train(Class ~ ., data = trainingSet,
              method = thresh_code,
              ## Minimize the distance to the perfect model
              metric = "Dist",
              maximize = FALSE,
              tuneLength = 20,
              ntree = 500,
              trControl = trainControl(method = "repeatedcv",
                                       repeats = 5,
                                       classProbs = TRUE,
                                       summaryFunction = fourStats))

mod1
plot_1 <- ggplot(mod1) 
ggsave(filename="plot_1.png", plot=plot_1)
####################################

library(reshape2)
metrics <- mod1$results[, c(2, 4:6)]
metrics <- melt(metrics, id.vars = "threshold",
                variable.name = "Resampled",
                value.name = "Data")
colnames(metrics) <- c('threshold',"Resampled",'Data')
plot_2 <- ggplot(data = metrics, aes(x = threshold, y = Data, color = Resampled)) +
    geom_line() +
    ylab("") + xlab("Probability Cutoff") +
    theme(legend.position = "top")
ggsave(filename="plot_2.png", plot=plot_2)
