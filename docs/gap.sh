#!/usr/bin/bash

module load ceuadmin/R

cd ~/R
Rscript -e 'setwd("gap");devtools::document()'

export version=$(awk '/Version/{print $2}' ~/R/gap/DESCRIPTION)
export R_HOME=/usr/local/Cluster-Apps/ceuadmin/R/4.4.1-icelake/lib64/R

R CMD build --compact-vignettes=both --md5 --resave-data --log gap
rm gap/src/*.so gap/src/*.o
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check --as-cran gap_${version}.tar.gz
