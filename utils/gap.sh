#!/usr/bin/bash

module load gcc/6 geos-3.6.2-gcc-5.4.0-vejexvy pcre2-10.20-gcc-5.4.0-tcuhtrb texlive

cd ~/R
Rscript -e 'setwd("gap");devtools::document()'

export version=1.3-2
R CMD build --compact-vignettes=both --md5 --resave-data --log gap
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check --as-cran gap_${version}.tar.gz
rm gap/src/*.so gap/src/*.o
