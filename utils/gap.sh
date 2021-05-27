#!/usr/bin/bash

function rd()
{
  cd gap
  Rscript -e "devtools::document()"
}

export version=1.2.3-2

R CMD build --compact-vignettes=both --md5 --resave-data gap
R CMD check gap_${version}.tar.gz
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
