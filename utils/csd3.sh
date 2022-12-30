#!/usr/bin/bash

module load gcc/6 geos-3.6.2-gcc-5.4.0-vejexvy pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
export gcc6=/usr/local/software/master/gcc/6
export version=4.2.2
export prefix=/rds-d4/user/$USER/hpc-work
export HPC_WORK=/rds/user/jhz22/hpc-work
cd ${prefix}
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz -
cd R-${version}
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${gcc6}/include:${HPC_WORK}/include LDFLAGS=-L${gcc6}/lib64:${gcc6}/lib:${HPC_WORK}/lib LIBS=-ltinfo
make
make install
ln -sf  ${prefix}/R-${version}/bin/R $HOME/bin/R
Rscript -e 'update.packages(checkBuilt=TRUE,ask=FALSE)'
