# R

[![pages-build-deployment](https://github.com/jinghuazhao/R/actions/workflows/pages/pages-build-deployment/badge.svg)](https://github.com/jinghuazhao/R/actions/workflows/pages/pages-build-deployment)

## A summary

Shown below are links to packages I have maintained[^packages], with individual files to be found from <https://github.com/jinghuazhao/R>.

See also documentation reformatted by datacamp, <https://www.rdocumentation.org/>: [gap](https://www.rdocumentation.org/packages/gap/),
 [gap.datasets](https://www.rdocumentation.org/packages/gap.datasets),
 [lmm](https://www.rdocumentation.org/packages/lmm),
 [pan](https://www.rdocumentation.org/packages/pan).
 [tdthap](https://www.rdocumentation.org/packages/tdthap).

**Packages**[^ctv] | [CRAN](http://cran.r-project.org) | Vignette | [GitHub](https://github.com/cran) | [R package documentation](https://rdrr.io/)
--------|---------------------------------------------|---------|---------------------|---------------------------------------------
**gaawr2** | <https://cran.r-project.org/package=gaawr2> | [gaawr2](https://jinghuazhao.github.io/R/vignettes/gaawr2.html) ([source](https://raw.githubusercontent.com/jinghuazhao/R/refs/heads/master/vignettes/gaawr2.Rmd), [R code](https://raw.githubusercontent.com/jinghuazhao/R/refs/heads/master/vignettes/gaawr2.R)), [web](https://jinghuazhao.github.io/R/vignettes/web.html) ([source](https://raw.githubusercontent.com/jinghuazhao/R/refs/heads/master/vignettes/web.Rmd), [R code](https://raw.githubusercontent.com/jinghuazhao/R/refs/heads/master/vignettes/web.R)) | <https://github.com/cran/gaawr2> | <https://rdrr.io/cran/gaawr2/>
**gap** | <https://cran.r-project.org/package=gap> | [gap](https://jinghuazhao.github.io/R/vignettes/gap.html)[^use] ([source](https://raw.githubusercontent.com/jinghuazhao/R/refs/heads/master/vignettes/gap.Rmd), [R code](https://raw.githubusercontent.com/jinghuazhao/R/refs/heads/master/vignettes/gap.R))   | <https://github.com/cran/gap> | <https://rdrr.io/cran/gap/>
 &nbsp; | &nbsp; | [jss](https://doi.org/10.18637/jss.v023.i08)
 &nbsp; | &nbsp; | [manual](https://jinghuazhao.github.io/R/vignettes/gap-manual.pdf)
 &nbsp; | &nbsp; | [shinygap](https://jinghuazhao.github.io/R/vignettes/shinygap.html)
**gap.datasets** | <https://cran.r-project.org/package=gap.datasets> | &nbsp; | <https://github.com/cran/gap.datasets> | <https://rdrr.io/cran/gap.datasets/>
**gap.examples** | &nbsp; | [gap.examples](https://jinghuazhao.github.io/R/vignettes/gap.examples.pdf)
 &nbsp;          | &nbsp; | [h2](https://jinghuazhao.github.io/R/vignettes/h2.pdf)
 &nbsp;          | &nbsp; | [pedtodot](https://jinghuazhao.github.io/R/vignettes/pedtodot.pdf)
 &nbsp;          | &nbsp; | [rnews](https://jinghuazhao.github.io/R/vignettes/rnews.pdf)
**lmm** | <https://cran.r-project.org/package=lmm> | [lmm](https://cran.r-project.org/web/packages/lmm/vignettes/lmm-tr.pdf) | <https://github.com/cran/lmm> | <https://rdrr.io/cran/lmm/>
**pan** | <https://cran.r-project.org/package=pan> | [pan](https://cran.r-project.org/web/packages/pan/vignettes/pan-tr.pdf) | <https://github.com/cran/pan> | <https://rdrr.io/cran/pan/>
**tdthap**  | <https://cran.r-project.org/package=tdthap> | [tdthap](https://cran.r-project.org/web/packages/tdthap/vignettes/tdthap-paper.pdf)| [https://github.com/cran/tdthap](https://github.com/cran/tdthap) | <https://rdrr.io/cran/tdthap/>
**kinship**[^kinship] | <https://cran.r-project.org/src/contrib/Archive/kinship/> | [kinship](https://jinghuazhao.github.io/R/vignettes/kinship.pdf) | <https://github.com/cran/kinship>

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

[^packages]: Contributed work

    > - Network Regression method in TWAS, [NeRiT](https://github.com/XiuyuanJin/NeRiT).
    > - `g.binread` function, [**GGIR**](https://cran.r-project.org/package=GGIR) & [GGIRread](https://cran.r-project.org/web/packages/GGIRread/index.html).
    > - [**ITHIM** injurymodel](https://github.com/ithim/injurymodel) through a hackathon at MRC.

[^ctv]: Featured CRAN task views:

    > - Genetics: **gap**, **tdthap**
    > - Meta-Analysis: **gap**
    > - Social science: **lmm**, **pan**

[^use]: Additional use can be found from a package vignette [pQTLtools](https://jinghuazhao.github.io/pQTLtools/articles/pQTLtools.html).

[^kinship]: Windows package [kinship_1.1.4.zip](kinship_1.1.4.zip) is built from [kinship_1.1.4.tar.gz](kinship_1.1.4.tar.gz) via <https://win-builder.r-project.org/>.
