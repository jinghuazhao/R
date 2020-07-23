L <- 53
INF <- Sys.getenv("INF")
r <- readLines(paste(INF,"doc","Olink.R",sep='/'),n=L)
write(r,file=paste(L))
source(paste(L))
save(hgTables,
     Cardiometabolic,
     Cell_Regulation,
     CVD_II,
     CVD_III,
     Development,
     Immune_Response,
     Immuno_Oncology,
     Inflammation,
     Metabolism,
     Neurology,
     Oncology_II,
     Organ_Damage,
     file='Olink_qPCR.rda',compress='xz')
unlink(paste(L))
unlink("inf1.csv")
unlink("inf2.csv")
