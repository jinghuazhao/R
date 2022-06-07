# R

## A summary

Links to packages I have maintained[^1] are as shown in the following table, while individual files can be found from GitHub.

**Packages**[^2] | [CRAN](http://cran.r-project.org) | [GitHub](https://github.com/cran) | [R package documentation](https://rdrr.io/)
--------|---------------------------------------------|------------------------------|---------------------------------------------
**gap** | [https://cran.r-project.org/package=gap](https://cran.r-project.org/package=gap)      | [https://github.com/cran/gap](https://github.com/cran/gap) | [https://rdrr.io/cran/gap/](https://rdrr.io/cran/gap/)
**gap.datasets** | [https://cran.r-project.org/package=gap.datasets](https://cran.r-project.org/package=gap.datasets) | [https://github.com/cran/gap.datasets](https://github.com/cran/gap.datasets) | [https://rdrr.io/cran/gap.datasets/](https://rdrr.io/cran/gap.datasets/)
**gap.examples** | 
**lmm** | [https://cran.r-project.org/package=lmm](https://cran.r-project.org/package=lmm)      | [https://github.com/cran/lmm](https://github.com/cran/lmm) | [https://rdrr.io/cran/lmm/](https://rdrr.io/cran/lmm/)
**pan** | [https://cran.r-project.org/package=pan](https://cran.r-project.org/package=pan)      | [https://github.com/cran/pan](https://github.com/cran/pan) | [https://rdrr.io/cran/pan/](https://rdrr.io/cran/pan/)
**tdthap**  | [https://cran.r-project.org/package=tdthap](https://cran.r-project.org/package=tdthap) | [https://github.com/cran/tdthap](https://github.com/cran/tdthap) | [https://rdrr.io/cran/tdthap/](https://rdrr.io/cran/tdthap/)
**kinship**[^3] | [https://cran.r-project.org/src/contrib/Archive/kinship/](https://cran.r-project.org/src/contrib/Archive/kinship/) | [https://github.com/cran/kinship](https://github.com/cran/kinship)

It contains packages **CGR** which is not available from CRAN and **kinship** which has additional update. 

## Vignettes

1. [gap](https://jinghuazhao.github.io/R/vignettes/gap.html). Package vignette[^4].
   * [jss](https://jinghuazhao.github.io/R/vignettes/jss.pdf). The JSS paper.
   * [shinygap](https://jinghuazhao.github.io/R/vignettes/shinygap.html). Shiny app on study designs.
2. [gap.examples](https://jinghuazhao.github.io/R/vignettes/gap.examples.pdf)
   * [h2](https://jinghuazhao.github.io/R/vignettes/h2.pdf). Polygenic modeling now contained in **gap.examples**.
   * [pedtodot](https://jinghuazhao.github.io/R/vignettes/pedtodot.pdf). Example usages of *pedtodot*.
   * [rnews](https://jinghuazhao.github.io/R/vignettes/rnews.pdf). *R News* application note.
3. [kinship](https://jinghuazhao.github.io/R/vignettes/kinship.pdf). Package vignette.
4. [lmm](https://cran.r-project.org/web/packages/lmm/vignettes/lmm-tr.pdf). Package vignette.
5. [pan](https://cran.r-project.org/web/packages/pan/vignettes/pan-tr.pdf). Package vignette.
6. [tdthap](https://cran.r-project.org/web/packages/tdthap/vignettes/tdthap-paper.pdf). Package vignette.

## Installation

You can install these packages either from CRAN, e.g.,
```r
install.packages("pan", repos="https://cran.r-project.org")
```
or GitHub, 
```r
library(devtools)
install_github("jinghuazhao/R/gap", build_vignettes = TRUE)
```

---

[^1]: Contributed work

    > - Network Regression method in TWAS, [NeRiT](https://github.com/XiuyuanJin/NeRiT).
    > - `g.binread` function in [**GGIR**](https://cran.r-project.org/package=GGIR) package.
    > - [**ITHIM** injurymodel](https://github.com/ithim/injurymodel) through a hackathon at MRC.

[^2]: Featured CRAN task views:

    > - Genetics: **gap**, **tdthap**
    > - Meta-Analysis: **gap**
    > - Social science: **lmm**, **pan**

[^3]: Windows package [kinship_1.1.4.zip](kinship_1.1.4.zip) is built from [kinship_1.1.4.tar.gz](kinship_1.1.4.tar.gz) via https://win-builder.r-project.org/.

[^4]: Additional use can be found from a package vignette [pQTLtools](https://jinghuazhao.github.io/pQTLtools/articles/pQTLtools.html).

