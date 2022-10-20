#!/usr/bin/bash

export wd=~/R
cd ${wd}/gap
Rscript -e "devtools::document()"
cd ${wd}/gap/inst/shinygap/

function rmd()
{
  Rscript -e "knitr::knit('README.Rmd')"
  pandoc README.md --citeproc --mathjax -s --self-contained -o index.html
}

cp README.Rmd ${wd}/gap/vignettes/shinygap.Rmd
cp shinygap.bib ${wd}/gap/vignettes/shinygap.bib
cd ${wd}

export version=1.3-1

R CMD build --compact-vignettes=both --md5 --resave-data --log gap
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check gap_${version}.tar.gz

rm ${wd}/gap/src/*.so ${wd}/gap/src/*.o
