# This is designed to find Machine Learning algorithm for predicting new compounds
# would be outliers or not.
# Class: 1 > outlier  0 > fitted 
# Input: 1. CountBB_0427.csv   (descirbe how SMILES is composed of building blocks) 
#        2. DFT_HOMO_BP86s_1percent_1outlier0fitted.csv (indicates which SMILES are outliers)

select_flavors <- "F18"






df_bb <- read.csv('BB_test_2.csv')
df_outlier <- read.csv('DFT_HOMO_BP86s_1percent_1outlier0fitted.csv')
list_1 <- df_bb$SMILES == df_outlier$SMILES
if(! FALSE %in% list_1){
    df_bb <- data.frame(df_bb,Class = df_outlier[,select_flavors])
}else{print("Two files have different SMILES!");stop()}
