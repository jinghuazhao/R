---
title: Shiny for Genetic Analysis Package (gap) Designs
subtitle: "Jing Hua Zhao"
output:
  bookdown::html_document2:
    toc: true
    toc_float: true
fontsize: 11pt
lot: true
lof: true
author:
  - Department of Public Health and Primary Care
  - University of Cambridge
  - Cambridge, UK
  - https://jinghuazhao.github.io/
date: '2023-03-02'
bibliography: '/rds/user/jhz22/hpc-work/R/gap/REFERENCES.bib'
csl: nature-genetics.csl
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Shiny for Genetic Analysis Package (gap) Designs}
  %\VignetteEncoding{UTF-8}
---

This is an initial attempt to enable easy calculation/visualization of study designs from R/gap which benchmarked relevant publications and eventually the app can produce more generic results.

One can run the app with R/gap installation as follows,

```r
setwd(file.path(find.package("gap"),"shinygap"))
library(shiny)
runApp()
```

Alternatively, one can run the app from source using `gap/inst/shinygap`. In fact, these are conveniently wrapped up as `runshinygap()` function.

To set the default parameters, some compromises need to be made, e.g., Kp=[1e-5, 0.4], MAF=[1e-3, 0.8], alpha=[1e-8, 0.05], beta=[0.01, 0.4]. The slider inputs provide upper bounds of parameters.

# Family-based study

This is a call to `fbsize()`.

# Population-based study

This is a call to `pbsize()`.

# Case-cohort study

This is a call to `ccsize()` whose `power` argument indcates power (TRUE) or sample size (FALSE) calculation.

# Two-stage case-control design

We implement it in function \texttt{tscc} whose format is
```
tscc(model, GRR, p1, n1, n2, M, alpha.genome, pi.samples, pi.markers, K)
```
which requires specification of disease model (multiplicative, additive, dominant, recessive), genotypic relative risk (GRR), the
estimated risk allele frequency in cases ($p_1$), total number of cases ($n_1$) total number of controls ($n_2$), total number of
markers ($M$), the false positive rate at genome level ($\alpha_\mathit{genome}$), the proportion of markers to be selected
($\pi_\mathit{markers}$, also used as the false positive rate at stage 1) and the population prevalence ($K$).

---

# Appendix: Theory

## A. Family-based and population-based designs

This is detailed in the package vignettes gap, <https://cran.r-project.org/web/packages/gap/vignettes/gap.html>, or jss, Zhao (2007).

## B. Case-cohort design

Our implemention is with respect to two aspects, Cai & Zeng (2004).

### B.1 Power

$$\Phi\left(Z_\alpha+\tilde{n}^\frac{1}{2}\theta\sqrt{\frac{p_1p_2p_D}{q+(1-q)p_D}}\right)$$

where $\alpha$ is the significance level, $\theta$ is the log-hazard ratio for
two groups, $p_j, j = 1, 2$, are the proportion of the two groups
in the population ($p_1 + p_2 = 1$), $\tilde{n}$ is the total number of subjects in the subcohort, $p_D$ is the proportion of the failures in
the full cohort, and $q$ is the sampling fraction of the subcohort.

### B.2 Sample size

$$\tilde{n}=\frac{nBp_D}{n-B(1-p_D)}$$ where $B=\frac{Z_{1-\alpha}+Z_\beta}{\theta^2p_1p_2p_D}$ and $n$ is the whole cohort size.

## C. Two-stage case-control design

Tests of allele frequency differences between cases and controls in a two-stage design are described here @skol06.
The usual test of proportions can be written as
$$z(p_1,p_2,n_1,n_2,\pi_{samples})=\frac{p_1-p_2}{\sqrt{\frac{p_1(1-p_1)}{2n_1\pi_{sample}}+\frac{p_2(1-p_2)}{2n_2\pi_{sample}}}}$$
where $p_1$ and $p_2$ are the allele frequencies, $n_1$ and $n_2$ are the sample sizes, $\pi_{samples}$ is the proportion of samples
to be genotyped at stage 1. The test statistics for stage 1, for stage 2 as replication and for stages 1 and 2 in a joint analysis
are then $z_1 = z(\hat p_1,\hat p_2,n_1,n_2,\pi_{samples})$, $z_2 = z(\hat p_1,\hat p_2,n_1,n_2,1-\pi_{samples})$,
$z_j = \sqrt{\pi_{samples}}z_1+\sqrt{1-\pi_{samples}}z_2$, respectively.
Let $C_1$, $C_2$, and $C_j$ be the thresholds for these statistics, the false positive rates can be obtained according to
$P(|z_1|>C_1)P(|z_2|>C_2,sign(z_1)=sign(z_2))$ and $P(|z_1|>C_1)P(|z_j|>C_j||z_1|>C_1)$ for replication-based and joint analyses, respectively.

# References
Zhao, J. H. [gap: Genetic analysis
package](https://doi.org/10.18637/jss.v023.i08).
*Journal of Statistical Software* **23**, 1--18 (2007).

Cai, J. & Zeng, D. [Sample size/power calculation for case-cohort
studies](https://doi.org/10.1111/j.0006-341X.2004.00257.x). *Biometrics*
**60**, 1015--24 (2004).

Skol, A. D., Scott, L. J., Abecasis, G. R. & Boehnke, M. [Joint analysis
is more efficient than replication-based analysis for two-stage
genome-wide association studies](https://doi.org/10.1038/ng1706). *Nat
Genet* **38**, 209--13 (2006).
