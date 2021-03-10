# R

This repository contains packages **CGR** which is not available from CRAN and **kinship** which has additional update. Links to packages I have maintained are as shown in the following table, with individual files listed by GitHub.

**Packages** | [CRAN](http://cran.r-project.org) | [GitHub](https://github.com/cran) | [R package documentation](https://rdrr.io/)
--------|---------------------------------------------|------------------------------|---------------------------------------------
**gap** | [https://cran.r-project.org/package=gap](https://cran.r-project.org/package=gap)      | [https://github.com/cran/gap](https://github.com/cran/gap) | [https://rdrr.io/cran/gap/](https://rdrr.io/cran/gap/)
**gap.datasets** | [https://cran.r-project.org/package=gap.datasets](https://cran.r-project.org/package=gap.datasets) | [https://github.com/cran/gap.datasets](https://github.com/cran/gap.datasets) | [https://rdrr.io/cran/gap.datasets/](https://rdrr.io/cran/gap.datasets/)
**gap.examples** | 
**lmm** | [https://cran.r-project.org/package=lmm](https://cran.r-project.org/package=lmm)      | [https://github.com/cran/lmm](https://github.com/cran/lmm) | [https://rdrr.io/cran/lmm/](https://rdrr.io/cran/lmm/)
**pan** | [https://cran.r-project.org/package=pan](https://cran.r-project.org/package=pan)      | [https://github.com/cran/pan](https://github.com/cran/pan) | [https://rdrr.io/cran/pan/](https://rdrr.io/cran/pan/)
**tdthap**  | [https://cran.r-project.org/package=tdthap](https://cran.r-project.org/package=tdthap) | [https://github.com/cran/tdthap](https://github.com/cran/tdthap) | [https://rdrr.io/cran/tdthap/](https://rdrr.io/cran/tdthap/)
**kinship**<sup>[1](#footnote1)</sup> | [https://cran.r-project.org/src/contrib/Archive/kinship/](https://cran.r-project.org/src/contrib/Archive/kinship/) | [https://github.com/cran/kinship](https://github.com/cran/kinship)

The [gap](https://jinghuazhao.github.io/R/vignettes/gap.pdf), [h2](https://jinghuazhao.github.io/R/vignettes/h2.pdf), [rnews](https://jinghuazhao.github.io/R/vignettes/rnews.pdf), [kinship](https://jinghuazhao.github.io/R/vignettes/kinship.pdf) in [vignettes](vignettes) and [pQTLtools](https://jinghuazhao.github.io/pQTLtools/).

Packages **gap** and **tdthap** are featured in [task view for genetics](https://cran.r-project.org/web/views/Genetics.html), while packages **lmm** and **pan** are featured in [task view for social sciences](https://cran.r-project.org/web/views/SocialSciences.html).

You can install these packages either from CRAN, e.g., 
```r
install.packages("pan", repos="https://cran.r-project.org")
```
or GitHub, 
```r
library(devtools)
install_github("jinghuazhao/R/gap", build_vignettes = TRUE)
```
I have implemented an earlier version of the `g.binread` function in [**GGIR**](https://cran.r-project.org/package=GGIR) package.

I have contributed to [**ITHIM** injurymodel](https://github.com/ithim/injurymodel) through a hackathon at MRC and recently to [cambridge-ceu](https://github.com/cambridge-ceu).

---

<a name="footnote1"><sup>1</sup></a> Windows package [kinship_1.1.4.zip](kinship_1.1.4.zip) is built from [kinship_1.1.4.tar.gz](kinship_1.1.4.tar.gz) via https://win-builder.r-project.org/.
