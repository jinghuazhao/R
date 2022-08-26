# R

## A summary

Shown below are links to packages I have maintained[^1], with individual files to be found from [https://github.com/jinghuazhao/R](https://github.com/jinghuazhao/R).

**Packages**[^2] | [CRAN](http://cran.r-project.org) | Vignette | [GitHub](https://github.com/cran) | [R package documentation](https://rdrr.io/)
--------|---------------------------------------------|---------|---------------------|---------------------------------------------
**gap** | [https://cran.r-project.org/package=gap](https://cran.r-project.org/package=gap) | [gap](https://jinghuazhao.github.io/R/vignettes/gap.html)[^3]   | [https://github.com/cran/gap](https://github.com/cran/gap) | [https://rdrr.io/cran/gap/](https://rdrr.io/cran/gap/)
 &nbsp; | &nbsp; | [jss](https://jinghuazhao.github.io/R/vignettes/jss.pdf)
 &nbsp; | &nbsp; | [shinygap](https://jinghuazhao.github.io/R/vignettes/shinygap.html)
**gap.datasets** | [https://cran.r-project.org/package=gap.datasets](https://cran.r-project.org/package=gap.datasets) | &nbsp; | [https://github.com/cran/gap.datasets](https://github.com/cran/gap.datasets) | [https://rdrr.io/cran/gap.datasets/](https://rdrr.io/cran/gap.datasets/)
**gap.examples** | &nbsp; | [gap.examples](https://jinghuazhao.github.io/R/vignettes/gap.examples.pdf)
 &nbsp;          | &nbsp; | [h2](https://jinghuazhao.github.io/R/vignettes/h2.pdf)
 &nbsp;          | &nbsp; | [pedtodot](https://jinghuazhao.github.io/R/vignettes/pedtodot.pdf)
 &nbsp;          | &nbsp; | [rnews](https://jinghuazhao.github.io/R/vignettes/rnews.pdf)
**lmm** | [https://cran.r-project.org/package=lmm](https://cran.r-project.org/package=lmm) | [lmm](https://cran.r-project.org/web/packages/lmm/vignettes/lmm-tr.pdf) | [https://github.com/cran/lmm](https://github.com/cran/lmm) | [https://rdrr.io/cran/lmm/](https://rdrr.io/cran/lmm/)
**pan** | [https://cran.r-project.org/package=pan](https://cran.r-project.org/package=pan) | [pan](https://cran.r-project.org/web/packages/pan/vignettes/pan-tr.pdf) | [https://github.com/cran/pan](https://github.com/cran/pan) | [https://rdrr.io/cran/pan/](https://rdrr.io/cran/pan/)
**tdthap**  | [https://cran.r-project.org/package=tdthap](https://cran.r-project.org/package=tdthap) | [tdthap](https://cran.r-project.org/web/packages/tdthap/vignettes/tdthap-paper.pdf)| [https://github.com/cran/tdthap](https://github.com/cran/tdthap) | [https://rdrr.io/cran/tdthap/](https://rdrr.io/cran/tdthap/)
**kinship**[^4] | [https://cran.r-project.org/src/contrib/Archive/kinship/](https://cran.r-project.org/src/contrib/Archive/kinship/) | [kinship](https://jinghuazhao.github.io/R/vignettes/kinship.pdf) | [https://github.com/cran/kinship](https://github.com/cran/kinship)

It contains packages **CGR** which is not available from CRAN and **kinship** which has additional update. 

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
    > - `g.binread` function, [**GGIR**](https://cran.r-project.org/package=GGIR) & [GGIRread](https://cran.r-project.org/web/packages/GGIRread/index.html).
    > - [**ITHIM** injurymodel](https://github.com/ithim/injurymodel) through a hackathon at MRC.

[^2]: Featured CRAN task views:

    > - Genetics: **gap**, **tdthap**
    > - Meta-Analysis: **gap**
    > - Social science: **lmm**, **pan**

[^3]: Additional use can be found from a package vignette [pQTLtools](https://jinghuazhao.github.io/pQTLtools/articles/pQTLtools.html).

[^4]: Windows package [kinship_1.1.4.zip](kinship_1.1.4.zip) is built from [kinship_1.1.4.tar.gz](kinship_1.1.4.tar.gz) via https://win-builder.r-project.org/.

