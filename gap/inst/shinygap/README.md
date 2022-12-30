---
title: Shiny for Genetic Analysis Package (gap) Designs
output:
  html_document:
    mathjax:  default
    fig_caption:  true
    toc: true
    section_numbering: true
bibliography: shinygap.bib
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

Alternatively, one can run the app from source using `gap/inst/shinygap`.

To set the default parameters, some compromises need to be made, e.g., Kp=[1e-5, 0.4], MAF=[1e-3, 0.8], alpha=[1e-8, 0.05], beta=[0.01, 0.4]. The slider inputs provide upper bounds of parameters.

## 1. Family-based study

This is a call to `fbsize()`.

## 2. Population-based study

This is a call to `pbsize()`.

## 3. Case-cohort study

This is a call to `ccsize()` whose `power` argument indcates power (TRUE) or sample size (FALSE) calculation.

## 4. Two-stage case-control design

This is a call to tscc().

---

# Appendix: Theory

## A. Family-based and population-based designs

See the R/gap package vignette jss or @zhao07.

## B. Case-cohort design

### B.1 Power

Following @cai04, we have
$$\Phi\left(Z_\alpha+\tilde{n}^\frac{1}{2}\theta\sqrt{\frac{p_1p_2p_D}{q+(1-q)p_D}}\right)$$

where $\alpha$ is the significance level, $\theta$ is the log-hazard ratio for
two groups, $p_j, j = 1, 2$, are the proportion of the two groups
in the population ($p_1 + p_2 = 1$), $\tilde{n}$ is the total number of subjects in the subcohort, $p_D$ is the proportion of the failures in
the full cohort, and $q$ is the sampling fraction of the subcohort.

### B.2 Sample size

$$\tilde{n}=\frac{nBp_D}{n-B(1-p_D)}$$ where $B=\frac{Z_{1-\alpha}+Z_\beta}{\theta^2p_1p_2p_D}$ and $n$ is the whole cohort size.

## C. Two-stage case-control design

In the notation of @skol06,

$$P(|z_1|>C_1)P(|z_2|>C_2,sign(z_1)=sign(z_2))$$ and $$P(|z_1|>C_1)P(|z_j|>C_j\,|\,|z_1|>C_1)$$
for replication-based and joint analyses, respectively; where $C_1$, $C_2$, and $C_j$
are threshoulds at stages 1, 2 replication and joint analysis,
$z_1 = z(p_1,p_2,n_1,n_2,\pi_{samples})$, $\,$
$z_2 = z(p_1,p_2,n_1,n_2,1-\pi_{samples})$, $\,$
$z_j = z_1 \sqrt{\pi_{samples}}+z_2\sqrt{1-\pi_{samples}}$.

# References

Cai, J., and D. Zeng. 2004. <span>“<a href="https://www.ncbi.nlm.nih.gov/15606422">Sample Size/Power Calculation for Case–Cohort Studies</a>.”</span> Journal Article. <em>Biometrics</em> 60 (4): 1015–24.

Skol, A. D., L. J. Scott, G. R. Abecasis, and M. Boehnke. 2006. <span>“Joint Analysis Is More Efficient Than Replication-Based Analysis for Two-Stage Genome-Wide Association Studies.”</span> <em>Nat Genet</em> 38 (2): 209–13.

Zhao, J. H. 2007. <span>“gap: Genetic Analysis Package.”</span><em>Journal of Statistical Software</em> 23 (8): 1–18. <a href="https://doi.org/10.18637/jss.v023.i08">https://doi.org/10.18637/jss.v023.i08</a>.
