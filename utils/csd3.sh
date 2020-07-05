#!/usr/bin/bash

module load gcc/6
module load pcre/8.38
module load texlive
wget https://cran.r-project.org/src/base/R-4/R-4.0.2.tar.gz
tar xvfz R-4.0.2.tar.gz
cd R-4.0.2
export prefix=/rds-d4/user/$USER/hpc-work
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${prefix}/include LDFLAGS=-L${prefix}/lib
make
make install

