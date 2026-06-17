set.seed(0)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.height = 8,
  fig.width = 8,
  fig.path = "gaawr2/",
  collapse = TRUE,
  comment = "#>",
  dev = "CairoPNG")

pkgs <- c("EnsDb.Hsapiens.v75","ensembldb","GMMAT","HardyWeinberg","MCMCglmm","SNPassoc","biomaRt",
          "gap","gap.datasets","haplo.stats","powerEQTL","R2jags","regress",
          "dplyr","ggplot2","httr","jsonlite","kableExtra","knitr","tidyr")
has_pkg <- function(x) requireNamespace(x, quietly = TRUE)
installed <- pkgs[vapply(pkgs, has_pkg, logical(1))]
missing <- setdiff(pkgs, installed)
if (length(missing)) {
  message("Missing packages: ", paste(missing, collapse = ", "))
}
invisible(lapply(installed, function(p)
  suppressMessages(library(p, character.only = TRUE))
))
sys_options <- options()
new_options <- options(digits=2)

desc <- read.dcf(system.file("DESCRIPTION", package="gaawr2"))
description <- desc[, "Description"]
description <- gsub(
  "(doi:[0-9\\.a-zA-Z]+/[0-9A-Z]+)",
  "[\\1](https://doi.org/\\1)",
  description
)
description <- gsub("\n+", " ", description)
knitr::asis_output(description)

print("Hello, world!\n")

## export message="Hello, world!"
## echo "print('$message')" > hello.R
## R CMD BATCH hello.R
## R --no-save -q < hello.R
## R --no-save -q <<END
## message <- Sys.getenv("message"); print(message)
## source("hello.R")
## END
## echo ${message} | \
## Rscript -e '
## message <- scan("stdin", what="", sep="\n", quiet=TRUE);
## write.table(message, col.names=FALSE, row.names=FALSE,
##             quote=FALSE)
## ' | \
## cat
## rm hello.*

class(iris)
dim(iris)
str(iris)
head(iris,1)
tail(iris,1)

options(new_options)
data(diabetes,package="gaawr2")

mean_values <- diabetes %>%
  dplyr::filter(CLASS %in% c("Y", "N", "P")) %>%
  dplyr::mutate(
    Gender = dplyr::recode(Gender, "F" = "Female", "M" = "Male"),
    CLASS = dplyr::recode(CLASS, "Y" = "Yes", "N" = "No", "P" = "Predicted")
  ) %>%
  dplyr::group_by(CLASS, Gender) %>%
  dplyr::select(AGE:BMI) %>%
  dplyr::summarize(dplyr::across(dplyr::everything(), \(x) mean(x, na.rm = TRUE)))
kableExtra::kbl(mean_values,caption="Mean value by gender and diabetes category") %>%
kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

mean_values_long <- mean_values %>%
  tidyr::pivot_longer(
    cols = AGE:BMI,
    names_to = "Variable",
    values_to = "Mean_Value"
  )
ggplot2::ggplot(mean_values_long,
  ggplot2::aes(x = Variable, y = Mean_Value, fill = CLASS)) +
  ggplot2::geom_col(position = ggplot2::position_dodge()) +
  ggplot2::facet_wrap(~ Gender) +
  ggplot2::labs(
    title = "",
    x = "Variable",
    y = "Mean Value",
    fill = "Diabetes Status" # Modified legend title for clarity
  ) +
  ggplot2::theme_bw() +
  ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
options(sys_options)

# MN blood group
SNP <- c(MM = 298, MN = 489, NN = 213)
HardyWeinberg::maf(SNP)
HardyWeinberg::HWTernaryPlot(SNP,region=0,grid=TRUE,markercol="blue")
HardyWeinberg::HWChisq(SNP, cc = 0, verbose = TRUE)
# Chromosome X
xSNP <- c(A=10, B=20, AA=30, AB=20, BB=10)
HardyWeinberg::HWChisq(xSNP,cc=0,x.linked=TRUE,verbose=TRUE)
# HLA/DQR
DQR <- gap.datasets::hla[,3:4]
a1 <- DQR[1]
a2 <- DQR[2]
GenotypeCounts <- HardyWeinberg::AllelesToTriangular(a1,a2)
kableExtra::kbl(GenotypeCounts,caption="Genotype distribution of DQR") %>%
kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
HardyWeinberg::HWPerm.mult(GenotypeCounts,nperm=300)
HardyWeinberg::HWStr(hla[,3:4],test="permutation",nperm=300)

data(asthma, package = "SNPassoc")
str(asthma, list.len=8)
knitr::kable(asthma[1:3,1:8],caption="First three records & two SNPs")
snpCols <- colnames(asthma)[6+(1:2)]
snps <- SNPassoc::setupSNP(data=asthma[snpCols], colSNPs=1:length(snpCols), sep="")
head(snps)
summary(snps, print=FALSE)
lapply(snps, head)
lapply(snps, summary)
SNPassoc::tableHWE(snps)

asthma.snps <- asthma %>%
               dplyr::rename(cc=casecontrol) %>%
               SNPassoc::setupSNP(colSNPs=(6+1):ncol(.), sep="")
# Model 1: Simple SNP association with BMI
SNPassoc::association(bmi ~ rs4490198, data = asthma.snps)

# Model 2: SNP association with case-control status
SNPassoc::association(cc ~ rs4490198, data = asthma.snps)

# Model 3: SNP association with covariates (country and smoke)
SNPassoc::association(cc ~ rs4490198 + country + smoke, data = asthma.snps)

# Model 4: SNP association with stratification by gender
SNPassoc::association(cc ~ rs4490198 + survival::strata(gender), data = asthma.snps)

# Model 5: SNP association with subset (only Spain)
SNPassoc::association(cc ~ rs4490198, data = asthma.snps, subset = country == "Spain")

# Model 6: Interaction between SNP (dominant model) and smoking
SNPassoc::association(cc ~ SNPassoc::dominant(rs4490198) * factor(smoke), data = asthma.snps)

# Model 7: Interaction between two SNPs (dominant model for rs4490198)
SNPassoc::association(cc ~ rs4490198 * factor(rs11123242), data = asthma.snps, model.interaction = "dominant")

# Association with the first three SNPs
snpsH <- names(asthma.snps)[6+(1:3)]
genoH <- SNPassoc::make.geno(asthma.snps, snpsH)
em <- haplo.stats::haplo.em(genoH, locus.label = snpsH, miss.val = c(0, NA))
haplo_table <- with(em,cbind(haplotype,hap.prob))
knitr::kable(haplo_table,caption="Haplotypes of the first three SNPs")
modH <- haplo.stats::haplo.glm(cc ~ genoH, data=asthma.snps,
                               family="binomial",
                               locus.label=snpsH,
                               allele.lev=attributes(genoH)$unique.alleles,
                               control = haplo.stats::haplo.glm.control(haplo.freq.min=0.05))
modH
SNPassoc::intervals(modH)

# Model comparison with / without haplotypes
mod.adj.ref <- glm(cc ~ smoke, data=asthma.snps, family="binomial")
mod.adj <- haplo.glm(cc ~ genoH + smoke, data=asthma.snps,
                 family="binomial",
                 locus.label=snpsH,
                 allele.lev=attributes(genoH3)$unique.alleles,
                 control = haplo.stats::haplo.glm.control(haplo.freq.min=0.05))
mod.adj
lrt.adj <- mod.adj.ref$deviance - mod.adj$deviance
pchisq(lrt.adj, mod.adj$lrt$df, lower=FALSE)

# Four variable slide windows over nine SNPs
snpsH <- labels(asthma.snps)[6+(1:9)]
genoH <- SNPassoc::make.geno(asthma.snps, snpsH)
haploH <- list()
for (i in 1:4) haploH[[i]] <- haplo.stats::haplo.score.slide(asthma.snps$cc, genoH,
                              trait.type="binomial",
                              n.slide=i,
                              locus.label=snpsH,
                              simulate=TRUE,
                              sim.control=haplo.stats::score.sim.control(min.sim=50,max.sim=100))

data(example,package="GMMAT")
attach(example)
model0 <- GMMAT::glmmkin(disease ~ age + sex, data = pheno, kins = GRM,
                         id = "id", family = binomial(link = "logit"))
model1 <- GMMAT::glmmkin(fixed = trait ~ age + sex, data = pheno, kins = GRM,
                         id = "id", family = gaussian(link = "identity"))
model2 <- GMMAT::glmmkin(fixed = trait ~ age + sex, data = pheno, kins = GRM,
                         id = "id", groups = "disease",
                         family = gaussian(link = "identity"))
snps <- c("SNP10", "SNP25", "SNP1", "SNP0")
geno.file <- system.file("extdata", "geno.bgen", package = "GMMAT")
samplefile <- system.file("extdata", "geno.sample", package = "GMMAT")
outfile <- "glmm.score.txt"
GMMAT::glmm.score(model0, infile = geno.file, BGEN.samplefile = samplefile,
                  outfile = outfile)
read.delim(outfile) |>
     head(n=4) |>
     knitr::kable(caption="Score tests under GLMM on four SNPs",digits=2)
unlink(outfile)
bed.file <- system.file("extdata", "geno.bed", package = "GMMAT") |>
            tools::file_path_sans_ext()
model.wald <- GMMAT::glmm.wald(fixed = disease ~ age + sex, data = pheno,
                               kins = GRM, id = "id", family = binomial(link = "logit"),
                               infile = bed.file, snps = snps)
knitr::kable(model.wald,caption="Wald tests under GLMM on four SNPs")
detach(example)

set.seed(1234567)
meyer <- within(gap.datasets::meyer,{
         y[is.na(y)] <- rnorm(length(y[is.na(y)]),mean(y,na.rm=TRUE),sd(y,na.rm=TRUE))
         g1 <- ifelse(generation==1,1,0)
         g2 <- ifelse(generation==2,1,0)
         id <- animal
         animal <- ifelse(!is.na(animal),animal,0)
         dam <- ifelse(!is.na(dam),dam,0)
         sire <- ifelse(!is.na(sire),sire,0)
     })
G <- gap::kin.morgan(meyer)$kin.matrix*2
r <- regress::regress(y~-1+g1+g2,~G,data=meyer)
r
with(r,gap::h2G(sigma,sigma.cov))
eps <- 0.001
y <- with(meyer,y)
x <- with(meyer,cbind(g1,g2))
ex <- gap::h2.jags(y,x,G,sigma.p=0.03,sigma.r=0.014,n.chains=1,n.iter=80)
kableExtra::kbl(ex$BUGSoutput$summary,digits=2,caption="MCMC results for the Meyer data") %>%
kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))

n.designs <- 6
designs <- 1:n.designs
N <- 50 * designs
n.grids <- 100
index <- 1:n.grids
grids <- index / n.grids
MAF <- seq(0.005, n.grids/2, by=0.5) / n.grids
plot(MAF,grids,type="n",ylab="Power")
mtext(expression(paste("(",alpha," = 0.05)")),1,line=4.5)
colors <- grDevices::hcl.colors(n.designs)
for (design in designs)
{
  power.SLR <- rep(NA,n.grids)
  for (j in index) power.SLR[j] <- powerEQTL::powerEQTL.SLR(MAF = MAF[j], FWER = 0.05, nTests = 240, slope = 0.13,
                                                            n = N[design], sigma.y = 0.13)
  lines(MAF,power.SLR,col=colors[design])
}
legend("bottomright", inset=.02, title="Sample size (N)", paste(N), col=colors, horiz=FALSE, cex=0.8, lty=designs)

data("EnsDb.Hsapiens.v75", package="EnsDb.Hsapiens.v75")
ensembldb::metadata(EnsDb.Hsapiens.v75)
genes <- ensembldb::genes(EnsDb.Hsapiens.v75)
head(genes)
transcripts_data <- ensembldb::transcripts(EnsDb.Hsapiens.v75)
head(transcripts_data)

gene_id <- "ENSG00000164308"
query_string = "
  query target($ensemblId: String!){
    target(ensemblId: $ensemblId){
      id
      approvedSymbol
      biotype
      geneticConstraint {
        constraintType
        exp
        obs
        score
        oe
        oeLower
        oeUpper
      }
      tractability {
        label
        modality
        value
      }
    }
  }
"
base_url <- "https://api.platform.opentargets.org/api/v4/graphql"
variables <- list("ensemblId" = gene_id)
post_body <- list(query = query_string, variables = variables)
r <- httr::POST(url=base_url, body=post_body, encode='json')
data <- iconv(r, "", "ASCII")
content <- jsonlite::fromJSON(data)
target <- content$data$target
scalar_fields <- data.frame(
  Field = c("ID", "Approved Symbol", "Biotype"),
  Value = c(target$id, target$approvedSymbol, target$biotype)
)
tractability_data <- target$tractability
kableExtra::kbl(scalar_fields,caption="(a) Basic Information") %>%
kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
kableExtra::kbl(target$geneticConstraint, caption="(b) Genetic Constraint Metrics") %>%
kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
kableExtra::kbl(tractability_data,caption="(c) Tractability Information") %>%
kableExtra::kable_styling(bootstrap_options = c("striped", "hover"), full_width = FALSE)

suggests <- read.dcf(file = system.file("DESCRIPTION", package = "gaawr2"), fields = c("Suggests"))
write.dcf(suggests)
