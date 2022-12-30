#!/usr/bin/bash

module load gcc/6 texlive
module load pcre2-10.20-gcc-5.4.0-tcuhtrb
module load geos-3.6.2-gcc-5.4.0-vejexvy

export wd=~/R
cd ${wd}/gap
Rscript -e "devtools::document()"
cd ${wd}/gap/inst/shinygap/
cp README.Rmd ${wd}/gap/vignettes/shinygap.Rmd
cp shinygap.bib ${wd}/gap/vignettes/shinygap.bib

# Rscript -e "knitr::knit('README.Rmd')"
# pandoc README.md --citeproc --mathjax -s --self-contained -o index.html

cd ${wd}
export version=1.3-2

R CMD build --compact-vignettes=both --md5 --resave-data --log gap
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check --as-cran gap_${version}.tar.gz

rm ${wd}/gap/src/*.so ${wd}/gap/src/*.o
