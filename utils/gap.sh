#!/usr/bin/bash

export wd=~/R
cd ${wd}/gap
Rscript -e "devtools::document()"
cd ${wd}/gap/inst/shinygap/
Rscript -e "knitr::knit('README.Rmd')"
pandoc README.md --citeproc --mathjax -s --self-contained -o index.html
cp index.html ${wd}/vignettes/shinygap.html
cp README.Rmd ${wd}/gap/vignettes/shinygap.Rmd
cp shinygap.bib ${wd}/gap/vignettes/shinygap.bib
cd ${wd}

export version=1.2.3-2

R CMD build --compact-vignettes=both --md5 --resave-data gap
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check gap_${version}.tar.gz
