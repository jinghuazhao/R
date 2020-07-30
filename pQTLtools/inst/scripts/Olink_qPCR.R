Lines <- 36
INF <- Sys.getenv("INF")
r <- readLines(paste(INF,"doc","Olink.R",sep='/'),n=Lines)
write(r,file=paste(Lines))
source(paste(Lines))
hgTables <- subset(hgTables,!grepl("hap",X.chrom)&!grepl("Un",X.chrom)&!grepl("random",X.chrom)&!grepl(";",geneName)&geneName!="")
save(hgTables,file="hgTables.rda",compress='xz')
save(Cardiometabolic,
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
unlink(paste(Lines))
unlink("inf1.csv")
unlink("inf2.csv")

combine <- function()
# Inflammation panel somehow has different header
{
  cols <- c(2,3,5)
  col_names <- c("UniProt","pctl10","pctl90")
  names(Cardiometabolic)[cols] <- col_names
  names(Cell_Regulation)[cols] <- col_names
  names(CVD_II)[cols] <- col_names
  names(CVD_III)[cols] <- col_names
  names(Development)[cols] <- col_names
  names(Immune_Response)[cols] <- col_names
  names(Immuno_Oncology)[cols] <- col_names
  names(Metabolism)[cols] <- col_names
  names(Neurology)[cols] <- col_names
  names(Oncology_II)[cols] <- col_names
  names(Organ_Damage)[cols] <- col_names
  Olink_qPCR <- rbind(
    data.frame(Panel="Cardiometabolic",Cardiometabolic),
    data.frame(Panel="Cell_Regulation",Cell_Regulation),
    data.frame(Panel="CVD_II",CVD_II),
    data.frame(Panel="CVD_III",CVD_III),
    data.frame(Panel="Development",Development),
    data.frame(Panel="Immune_Response",Immune_Response),
    data.frame(Panel="Immuno_Oncology",Immuno_Oncology),
    data.frame(Panel="Metabolism",Metabolism),
    data.frame(Panel="Neurology",Neurology),
    data.frame(Panel="Oncology_II",Oncology_II),
    data.frame(Panel="Organ_Damage",Organ_Damage)
  )
  save(Olink_qPCR,Inflammation,file='Olink_qPCR.rda',compress='xz')
}
