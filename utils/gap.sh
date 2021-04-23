#!/usr/bin/bash

export version=1.2.3-2
R CMD build gap
R CMD check gap_${version}.tar.gz
R CMD INSTALL gap_${version}.tar.gz
