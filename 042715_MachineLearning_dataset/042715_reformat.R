
inputFile <- 'BB_test.csv'
temporaryFile <- 'BB_test_2.csv'


con  <- file(inputFile, open = "r")
con2 <- file(temporaryFile, open = "w")

head_con2 <- 'SMILES,B01,B02,B03,B04,B05,B06,B07,B08,B09,B10,B11,B12,B13,B14,B15,B16,B17,B18,B19,B20,B21,B22,B23,B24,B25,B26'
write(head_con2,file = con2)
while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    oneLine <- gsub("\"", "", oneLine)     #strip out the percent signs
    oneLine <- gsub("\\[", "", oneLine)   #strip out the dollar signs
    oneLine <- gsub("\\]", "", oneLine)   #strip out the dollar signs
    cat(oneLine, file = con2, sep = "\n") #spit the line back out
} 

close(con)
close(con2)

# SMILES_BB <- read.csv('BB_test_2.csv')