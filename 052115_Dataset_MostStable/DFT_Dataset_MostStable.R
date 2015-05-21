

#1 = BP86/SVP    2 = B3LYP/SVP  3 = PBE0/SVP 4 = BH&HLYP/SVP 5 = M06-2X/SVP 
#6 = HF/SVP      7 = BP86/TZVP   8 = B3LYP/TZVP 9 = PBE0/TZVP 
df_1 <- read.table("testlr.dat",skip = 2,strip.white = TRUE,sep = ",",
                   col.names=c("ID","SMILES","f1","f2","f3","f4","f5","f6","f7","f8","f9"),
                   na.strings=c("","NONE"))