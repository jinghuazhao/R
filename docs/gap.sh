#!/usr/bin/bash

module load gcc/9 geos-3.6.2-gcc-5.4.0-vejexvy pcre2-10.20-gcc-5.4.0-tcuhtrb texlive

cd ~/R
Rscript -e 'setwd("gap");devtools::document()'

export version=$(awk '/Version/{print $2}' ~/R/gap/DESCRIPTION)
R CMD build --compact-vignettes=both --md5 --resave-data --log gap
rm gap/src/*.so gap/src/*.o
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check --as-cran gap_${version}.tar.gz
