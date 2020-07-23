
HOME <- Sys.getenv("HOME")
SomaLogic160410 <- read.delim(paste(HOME,"SomaLogic","doc","SOMALOGIC_Master_Table_160410_1129info.tsv",sep="/"))
save(SomaLogic160410, file="SomaLogic160410.rda", compress='xz')
