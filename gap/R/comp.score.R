kin.coef <- function(n.parents, n.sibs)
{
    n.person <- n.parents+n.sibs
    kin.coef <- matrix(0.25, ncol=n.person, nrow=n.person)
    diag(kin.coef) <- 0.5
    if (n.parents == 2) kin.coef[1,2] <- kin.coef[2,1] <- 0
    kin.coef
}        

vec <- function(A)  A[upper.tri(A,diag=FALSE)]

score.core <- function(b1, b2, r)
# b1 and b2 are standardized to have variance 1
# r = corr (b1, b2)
{    
    S1 <- max(b1, 0)^2
    S2 <- NULL
    if ((b1 > r*b2) & (b2 >= r*b1)) S2 <- (b1^2-2*r*b1*b2+b2^2)/(1-r^2)
    else if ((b1 >= 0) & (b2 < r*b1)) S2 <- b1^2
    else if ((b1 <= r*b2) & (b2 > 0)) S2 <- b2^2
    else S2 <- 0

    k <- acos(r)/(2*pi)
    p1 <- 0.5-0.5*pchisq(S1, df = 1)
    p2 <- 1-((0.5-k)+0.5*pchisq(S2, df=1)+k*pchisq(S2,df=2))
    c(round(S1,3), round(S2,3), round(p1,4), round(p2,4))
}

reorder <- function(familyID, pair.id)
# reorder the sib ID so the one with smaller ID comes first
# each sib ID is prefixed by its family ID
# familyID can be either a scaler or a vector
{
    pairs <- do.call("rbind", strsplit(as.character(pair.id), ","))
    p1 <- as.numeric(pairs[,1])
    p2 <- as.numeric(pairs[,2])
    ttt1 <- paste(familyID, pmin(p1, p2), sep=".")
    ttt2 <- paste(familyID, pmax(p1, p2), sep=".")
    paste(ttt1, ttt2, sep=",")
}

pair.data <- function(pheno,mean,h2,var)
# return vector v
{
    # get rid of person with missing trait value
    pheno <- pheno[!is.na(pheno$phenotype),] 
    # put parants ahead of sibs
    pheno <- pheno[order(pheno$family, pheno$father, pheno$mother),] 

    pairID <- allv <- NULL
    for (family in unique(pheno$family))
    {
        target <- (pheno$family == family)
        ttt <- pheno[target,]
        # a person is either a founder or a non-founder
        n.parents <- sum(ttt$father == 0 & ttt$mother == 0) 
        n.sibs <- sum(ttt$father!=0 | ttt$mother !=0)
        if (n.sibs < 2) next
        
        # no shared components 
        Sigma0 <- 2*h2*kin.coef(n.parents,n.sibs)+(1-h2)*diag(n.parents+n.sibs)    
        Sigma0 <- var*Sigma0
        Sigma0.1 <- solve(Sigma0)
        
        w <- as.vector(Sigma0.1 %*% (ttt$phenotype-mean))
        Wanted <- (n.parents+1):(n.parents+n.sibs)
        v <- vec((outer(w,w)-Sigma0.1)[Wanted,Wanted])
        names <- as.character(ttt$person)
        pair.id <- vec(outer(names,names, FUN=paste, sep=",")[Wanted,Wanted])

        n.pairs <- length(v)
        pairID <- c(pairID, reorder(family,pair.id))
        allv <- c(allv, v) 
    }
    data.frame(pairID=pairID, v=allv)
}

score <- function(chrom.pos, pair.data, ibd)
{
        ibddata <- ibd[ibd$pos == chrom.pos & ibd$prior.Z0 == 0.25 &
                       ibd$prior.Z1 == 0.5 & ibd$prior.Z2 == 0.25,
                       c("pedigree", "pair", "Z1", "Z2")]
#        ibddata <- ibd[ibd$pos == chrom.pos, c("pedigree","pair","Z1","Z2")]
        newdata <- data.frame(pairID=reorder(ibddata$pedigree, ibddata$pair),
                   epi=ibddata$Z1/2+ibddata$Z2, pi2=ibddata$Z2)
        joint.data <- merge(pair.data, newdata, by="pairID", all=FALSE)

        b1 <- sum((joint.data$epi - 0.5) * joint.data$v)
        b2 <- sum((joint.data$pi2 - 0.25) * joint.data$v)
        Vepi <- var(cbind(joint.data$epi, joint.data$pi2))
        Vb <- Vepi*sum(joint.data$v^2)
        b1 <- b1/sqrt(Vb[1,1])
        b2 <- b2/sqrt(Vb[2,2])
        r <- Vb[1,2]/sqrt(Vb[1,1]*Vb[2,2])

        c(chrom.pos, score.core(b1, b2, r))
}
################################### End of the main part ###################################


############################ Function for end users #################################

#' score statistics for testing genetic linkage of quantitative trait
#'
#' The function empirically estimate the variance of the score functions.
#' The variance-covariance matrix consists of two parts: the additive
#' part and the part for the individual-specific environmental effect.
#' Other reasonable decompositions are possible. 
#'
#' This program has the following improvement over "score.r":
#'
#'   1. It works with selected nuclear families
#'
#'   2. Trait data on parents (one parent or two parents), if available,
#'      are utilized.
#'
#'   3. Besides a statistic assuming no locus-specific dominance effect,
#'      it also computes a statistic that allows for such effect.
#'      It computes two statistics instead of one.  
#'
#' Function "merge" is used to merge the IBD data for a pair with the
#' transformed trait data (i.e., \eqn{w_kw_l}).
#'
#' @param ibddata The output file from GENEHUNTER using command "dump ibd". The default file name is \eqn{ibd_dist.out}.
#' @param phenotype The file of pedigree structure and trait value.  
#' The default file name is "pheno.dat". Columns (no headings) are: 
#' family ID, person ID, father ID, mother ID, gender, trait value,
#' where Family ID and person ID must be numbers, not characters.
#' Use character "NA" for missing phenotypes.
#' @param mean (population) mean of the trait, with a default value of 0.
#' @param var (population) variance of the trait, with a default value of 1.
#' @param h2 heritability of the trait, with a default value of 0.3.
#'
#' @export
#' @return a matrix with each row containing the location and the statistics and their p-values.
#'
#' @references
#' Kruglyak L, Daly MJ, Reeve-Daly MP, Lander ES (1996) Parametric and Nonparametric 
#' linkage analysis: a unified multipoint approach. Am J Hum Genet 58:1347-1363
#'
#' Kruglyak L, Lander ES (1998) Faster multipoint linkage analysis using Fourier transforms
#' J Comp Bio 1998 5:1-7
#'
#' Wang K (2005) A likelihood approach for quantitative-trait-locus mapping with
#' selected pedigrees. Biometrics 61:465-473
#'
#' @examples
#' \dontrun{
#' # An example based on GENEHUNTER version 2.1, with quantitative trait data in file
#' # "pheno.dat" generated from the  standard normal distribution. The following
#' # exmaple shows that it is possible to automatically call GENEHUNTER using R
#' # function "system".
#'
#' cwd <- getwd()
#' cs.dir <- file.path(find.package("gap"),"tests/comp.score")
#' setwd(cs.dir)
#' dir()
#' # system("gh < gh.inp")
#' cs.default <- comp.score()
#' setwd(cwd)
#' }
#'
#' @author Yingwei Peng, Kai Wang
#' @note Adapt from score2.r.
#' @keywords htest

comp.score <- function(ibddata="ibd_dist.out", phenotype="pheno.dat", mean=0, var=1, h2=0.3)
{
        ibd <- read.table(ibddata, skip = 1, col.names = c("pos", "pedigree",
               "pair", "prior.Z0", "prior.Z1", "prior.Z2", "Z0", "Z1", "Z2"),
               colClasses=c("numeric", "integer", "character", "numeric",
               "numeric", "numeric", "numeric", "numeric", "numeric"),
               comment.char = "")
        pheno <- read.table(phenotype, col.names = c("family",
                "person", "father", "mother", "gender", "phenotype"))
        paired <- pair.data(pheno,mean,h2,var)
        
        result <- NULL
        for (i in unique(ibd$pos)) result <- rbind(result, score(i, paired, ibd))
        dimnames(result) <- 
           list(NULL, c("pos", "stat.S1", "stat.S2", "p.value.S1", "p.value.S2"))
        result
}
