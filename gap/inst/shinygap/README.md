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
