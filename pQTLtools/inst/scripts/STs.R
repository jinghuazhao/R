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
st6.1 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(1:11:13), rows=c(3:167))
st6.2 <- openxlsx::read.xlsx(xlsx, sheet=6, colNames=TRUE, skipEmptyRows=TRUE,
                           cols=c(14:20), rows=c(4:167))
st6 <- cbind(st6.1,st6.2)
replicates <- merge(st4[,c(1:10,26:28)],st6[,c(1:10,17:20)],
                    by=c("Locus.ID","UniProt","Chr","Pos","SOMAmer.ID"))

# However, it is unclear UniProts in ST6 were selected from which of the panels
INF <- Sys.getenv("INF")
INF1_merge <- merge(inf1,
                    read.delim(file.path(INF,"work","INF1.merge-rsid"),as.is=TRUE),
                    by="prot")
INF1_uniprot <- unique(with(INF1_merge,uniprot))
options(width=250)
subset(replicates,UniProt %in% INF1_uniprot)
table(subset(replicates,UniProt %in% INF1_uniprot)$Replicates)

# side information on cvd2, cvd3, inf1
olink <- scan(file.path(INF,"doc","olink.prot.list.txt"),"")
olink_uniprot <- unlist(lapply(strsplit(olink,"___"),'[[',2))
dim(subset(replicates,UniProt %in% olink_uniprot))

