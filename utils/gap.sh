#!/usr/bin/bash

R CMD build gap
R CMD check gap_1.2.3.tar.gz
R CMD INSTALL gap_1.2.3.tar.gz
cp ~/hpc-work/R/gap/doc/gap.pdf vignettes/

