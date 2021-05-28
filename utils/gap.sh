#!/usr/bin/bash

function rd()
{
  cd gap
  Rscript -e "devtools::document()"
  cd -
  utils/st.sh
  Rscript -e "devtools::install_github(\"jinghuazhao/R/gap\",build_vignettes=TRUE,force=TRUE)"
}

export version=1.2.3-2

R CMD build --compact-vignettes=both --md5 --resave-data gap
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check gap_${version}.tar.gz
