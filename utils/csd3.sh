#!/usr/bin/bash

module load gcc/6
module load pcre/8.38
module load texlive
wget https://cran.r-project.org/src/base/R-4/R-4.0.3.tar.gz
tar xvfz R-4.0.3.tar.gz
cd R-4.0.3
export prefix=/rds-d4/user/$USER/hpc-work
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${prefix}/include LDFLAGS=-L${prefix}/lib
make
make install
cd $HOME/bin
ln -sf  /rds-d4/$HOME/hpc-work/R-4.0.3/bin/R
