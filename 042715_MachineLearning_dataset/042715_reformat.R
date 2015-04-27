
inputFile <- 'BB_test.csv'
temporaryFile <- 'BB_test_2.csv'


con  <- file(inputFile, open = "r")
con2 <- file(temporaryFile, open = "w")

head_con2 <- 'SMILES,B01,B02,B03,B04,B05,B06,B07,B08,B09,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B21,B22,B23,B24,B25,B26'
write(head_con2,file = con2)
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    position <- regexpr(",",oneLine) 
    str_1 <- substr(oneLine,start = 1L,stop = position-1)
    str_2 <- substr(oneLine,start = position,stop = 10000L)
    str_2 <- gsub("\"", "", str_2)     #strip out the  \" signs
    str_2 <- gsub("\\[", "", str_2)   #strip out the  \\[ signs
    str_2 <- gsub("\\]", "", str_2)   #strip out the  \\]signs
    oneLine <- paste(str_1,str_2,sep = "")
    cat(oneLine, file = con2, sep = "\n") #spit the line back out
} 

close(con)
close(con2)

SMILES_BB <- read.csv('BB_test_2.csv')