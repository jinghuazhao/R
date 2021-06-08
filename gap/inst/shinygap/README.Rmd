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

This is an initial attempt to enable easy calculation/visualization of three study designs.

One can run the app with R/gap installation as follows,

```r
setwd(file.path(find.package("gap"),"shinygap"))
library(shiny)
runApp()
```

Alternatively, one can run the app from source using `gap/inst/shinygap`.

## Family-based study

This is a call to `gap::fbsize()`. Note in particular that `alpha` parameter actually contains three elements and the 3rd of which is for ASP+TDT.

## Population-based study

This is a call to `gap::pbsize()`.

## Case-cohort study

This is a call to `gap::ccsize()`. Note that its `power` argument indcates power (TRUE) or sample size (FALSE) calculation.

To set the default parameters, some compromises need to be made, e.g., Kp=[1e-5, 0.4], MAF=[1e-3, 0.8], alpha=[5e-8, 0.05], beta=[0.01, 0.4]. The slider input provides upper bounds of a particular parameter.

It appears that `fbsize()` does not handle `gamma` well. Moreover, the `power` argument of `ccsize()` in gap 1.2.3-1 only accepts either NULL (default) or a fixed value `(1-beta)`..

---

# Appendix: Theory

## Family-based and population-based designs

See the R/gap package vignette jss or @zhao07.

## Case-cohort design

### Power

Following @cai04, we have
$$\Phi\left(Z_\alpha+\tilde{n}^\frac{1}{2}\theta\sqrt{\frac{p_1p_2p_D}{q+(1-q)p_D}}\right)$$

where $\alpha$ is the significance level, $\theta$ is the log-hazard ratio for
two groups, $p_j, j = 1, 2$, are the proportion of the two groups
in the population ($p_1 + p_2 = 1$), $\tilde{n}$ is the total number of subjects in the subcohort, $p_D$ is the proportion of the failures in
the full cohort, and $q$ is the sampling fraction of the subcohort.

### Sample size

$$\tilde{n}=\frac{nBp_D}{n-B(1-p_D)}$$ where $B=\frac{Z_{1-\alpha}+Z_\beta}{\theta^2p_1p_2p_D}$ and $n$ is the whole cohort size.

# References
