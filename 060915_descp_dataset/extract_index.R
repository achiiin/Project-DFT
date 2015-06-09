# df_1 <- read.csv('DFT_HOMO_BP86s_multi_1outlier0fitted_F18_V2.2.2.csv')
# df_2 <- read.table('uni_SMILES.smi')

df_1 <- read.csv('test_1.csv')
df_2 <- read.table('uni_SMILES_300.txt')
colnames(df_2) <- "SMILES"

df_2$No. <- 1:nrow(df_2)

df_3 <- merge(df_1,df_2,by="SMILES")
bad <- is.na(df_3$F18)
df_3 <- df_3[!bad,]

df_4 <- df_3[ , which(names(df_3) %in% c("No.","F18"))]
df_4 <- df_4[c("No.","F18")]
df_4 <- df_4[order(df_4[,"No."]),]
write.csv(df_4,file="F18_ID_class_v1.0.csv",row.names=FALSE,quote=FALSE)