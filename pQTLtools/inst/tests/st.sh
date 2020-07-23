#!/usr/bin/bash

R CMD build --resave-data pQTLtools
R CMD check --as-cran pQTLtools_0.1.tar.gz
R CMD INSTALL pQTLtools_0.1.tar.gz

# also work:
# R CMD INSTALL pQTLtools

