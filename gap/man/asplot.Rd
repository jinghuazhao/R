\name{asplot}
\alias{asplot}
\title{Regional association plot}
\usage{asplot(locus, map, genes, flanking=1e3, best.pval=NULL, sf=c(4,4),
              logpmax=10, pch=21)}
\description{
This function obtains regional association plot for a particular locus, based on
the information about recombinatino rates, linkage disequilibria between the
SNP of interest and neighbouring ones, and single-point association tests p values.

Note that the best p value is not necessarily within locus in the original design.
}

\arguments{
 \item{locus}{Data frame with columns c("CHR", "POS", "NAME", "PVAL", "RSQR") containing association results}
 \item{map}{Genetic map, i.e, c("POS","THETA","DIST")}
 \item{genes}{Gene annotation with columns c("START", "STOP", "STRAND", "GENE")}
 \item{flanking}{Flanking length}
 \item{best.pval}{Best p value for the locus of interest}
 \item{sf}{scale factors for p values and recombination rates, smaller values are necessary for gene dense regions}
 \item{logpmax}{Maximum value for -log10(p)}
 \item{pch}{Plotting character for the SNPs to be highlighted, e.g., 21 and 23 refer to circle and diamond}
}

\section{References}{
DGI. Whole-genome association analysis identifies novel loci for type 2 diabetes and triglyceride levels.
Science 2007;316(5829):1331-6
}

\examples{
\dontrun{
require(gap.datasets)
asplot(CDKNlocus, CDKNmap, CDKNgenes)
title("CDKN2A/CDKN2B Region")
asplot(CDKNlocus, CDKNmap, CDKNgenes, best.pval=5.4e-8, sf=c(3,6))

## NCBI2R

options(stringsAsFactors=FALSE)
p <- with(CDKNlocus,data.frame(SNP=NAME,PVAL))
hit <- subset(p,PVAL==min(PVAL,na.rm=TRUE))$SNP

library(NCBI2R)
# LD under build 36
chr_pos <- GetSNPInfo(with(p,SNP))[c("chr","chrpos")]
l <- with(chr_pos,min(as.numeric(chrpos),na.rm=TRUE))
u <- with(chr_pos,max(as.numeric(chrpos),na.rm=TRUE))
LD <- with(chr_pos,GetLDInfo(unique(chr),l,u))
# We have complaints; a possibility is to get around with 
# ftp://ftp.ncbi.nlm.nih.gov/hapmap/ld_data/2009-04_rel27/
# https://ftp.hapmap.org/hapmap/
hit_LD <- subset(LD,SNPA==hit)
hit_LD <- within(hit_LD,{RSQR=r2})
info <- GetSNPInfo(p$SNP)
haldane <- function(x) 0.5*(1-exp(-2*x))
locus <- with(info, data.frame(CHR=chr,POS=chrpos,NAME=marker,
                    DIST=(chrpos-min(chrpos))/1000000,
                    THETA=haldane((chrpos-min(chrpos))/100000000)))
locus <- merge.data.frame(locus,hit_LD,by.x="NAME",by.y="SNPB",all=TRUE)
locus <- merge.data.frame(locus,p,by.x="NAME",by.y="SNP",all=TRUE)
locus <- subset(locus,!is.na(POS))
ann <- AnnotateSNPList(p$SNP)
genes <- with(ann,data.frame(ID=locusID,CLASS=fxn_class,PATH=pathways,
                             START=GeneLowPoint,STOP=GeneHighPoint,
                             STRAND=ori,GENE=genesymbol,BUILD=build,CYTO=cyto))
attach(genes)
ugenes <- unique(GENE)
ustart <- as.vector(as.table(by(START,GENE,min))[ugenes])
ustop <- as.vector(as.table(by(STOP,GENE,max))[ugenes])
ustrand <- as.vector(as.table(by(as.character(STRAND),GENE,max))[ugenes])
detach(genes)
genes <- data.frame(START=ustart,STOP=ustop,STRAND=ustrand,GENE=ugenes)
genes <- subset(genes,START!=0)
rm(l,u,ugenes,ustart,ustop,ustrand)
# Assume we have the latest map as in CDKNmap
asplot(locus,CDKNmap,genes)
}
}
\author{Paul de Bakker, Jing Hua Zhao, Shengxu Li}
\keyword{hplot}
