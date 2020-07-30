single_panel <- function()
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

all_panels <- function()
# Inflammation panel somehow has different header
{
  cols <- c(1,2)
  col_names <- c("Target","UniProt")
  names(Cardiometabolic)[cols] <- col_names
  names(Cell_Regulation)[cols] <- col_names
  names(CVD_II)[cols] <- col_names
  names(CVD_III)[cols] <- col_names
  names(Development)[cols] <- col_names
  names(Immune_Response)[cols] <- col_names
  names(Immuno_Oncology)[cols] <- col_names
  names(Inflammation)[cols] <- col_names
  names(Metabolism)[cols] <- col_names
  names(Neurology)[cols] <- col_names
  names(Oncology_II)[cols] <- col_names
  names(Organ_Damage)[cols] <- col_names
  Olink_qPCR <- rbind(
    data.frame(Panel="Cardiometabolic",Cardiometabolic[cols]),
    data.frame(Panel="Cell_Regulation",Cell_Regulation[cols]),
    data.frame(Panel="CVD_II",CVD_II[cols]),
    data.frame(Panel="CVD_III",CVD_III[cols]),
    data.frame(Panel="Development",Development[cols]),
    data.frame(Panel="Immune_Response",Immune_Response[cols]),
    data.frame(Panel="Immuno_Oncology",Immuno_Oncology[cols]),
    data.frame(Panel="Inflammation",Inflammation[cols]),
    data.frame(Panel="Metabolism",Metabolism[cols]),
    data.frame(Panel="Neurology",Neurology[cols]),
    data.frame(Panel="Oncology_II",Oncology_II[cols]),
    data.frame(Panel="Organ_Damage",Organ_Damage[cols])
  )
}

Lines <- 36
INF <- Sys.getenv("INF")
r <- readLines(paste(INF,"doc","Olink.R",sep='/'),n=Lines)
write(r,file=paste(Lines))
source(paste(Lines))
hgTables <- subset(hgTables,!grepl("hap",X.chrom)&!grepl("Un",X.chrom)&!grepl("random",X.chrom)&!grepl(";",geneName)&geneName!="")
save(hgTables,file="hgTables.rda",compress='xz')
unlink(paste(Lines))
unlink("inf1.csv")
unlink("inf2.csv")

p12 <- all_panels()

library(biomaRt)
listEnsemblArchives()

# hg19/GRCh37
hg19 <- useMart(biomart= "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", host = "http://apr2020.archive.ensembl.org")
listAttributes(hg19)
hg19.bm <- getBM(attributes = c('uniprotswissprot', 'hgnc_symbol','chromosome_name', 'start_position', 'end_position'),
                 filters = 'uniprotswissprot',
                 values = unique(with(p12,UniProt)),
                 mart = hg19)
hg19.bed <- subset(hg19.bm,chromosome_name%in%c(1:22,'X','Y'))
names(hg19.bed) <- c("UniProt","gene","chr","start","end")
Olink_qPCR <- merge(p12,hg19.bed,by="UniProt",all=TRUE)
save(Olink_qPCR,file='Olink_qPCR.rda',compress='xz')
