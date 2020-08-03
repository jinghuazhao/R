# INTERVAL SomaLogic results
dir <- 'https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-018-0175-2/MediaObjects/'
file <- '41586_2018_175_MOESM4_ESM.xlsx'
xlsx <- paste0(dir,file)
st4.1 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                             cols=c(1:16,26:28), rows=c(5:1986))
st4.2 <- openxlsx::read.xlsx(xlsx, sheet=4, colNames=TRUE, skipEmptyRows=TRUE,
                             cols=c(17:25,29:31), rows=c(6:1986))
st4 <- cbind(st4.1,st4.2)
save(st4,file='st4.rda',compress='xz')
st5 <- openxlsx::read.xlsx(xlsx, sheet=5, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:19), rows=c(3:2746))
st6 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:20), rows=c(3:167))

