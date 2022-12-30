#!/usr/bin/bash

# Installation of R-devel
R-devel CMD INSTALL --configure-args="
 CC=\"/usr/bin/gcc\" \
 CXX=\"/usr/g++\" \
 FC=\"/usr/bin/gfortran\" \
 CFLAGS=\"-g -O2 -Wall -pedantic -mtune=native\" \
 FFLAGS=\"-g -O2 -mtune=native -Wall -pedantic\" \
 CXXFLAGS=\"-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses\" \
 LDFLAGS=\"-L/usr/lib64" $1

# GitHub installation example for R/gap
Rscript -e "devtools::install_github('jinghuazhao/R/gap',build_vignettes=TRUE,force=TRUE)"
cp ~/hpc-work/R/gap/doc/gap.pdf gap/vignettes/
