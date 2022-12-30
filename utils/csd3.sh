#!/usr/bin/bash

module load gcc/6 texlive
module load pcre2-10.20-gcc-5.4.0-tcuhtrb
module load geos-3.6.2-gcc-5.4.0-vejexvy

export version=4.2.2
export prefix=/rds-d4/user/$USER/hpc-work
export HPC_WORK=/rds/user/jhz22/hpc-work
cd ${prefix}
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz -
cd R-${version}
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${HPC_WORK}/include LDFLAGS=-L${HPC_WORK}/lib LIBS=-ltinfo
make
make install
cd $HOME/bin
ln -sf  ${prefix}/R-${version}/bin/R
Rscript -e 'update.packages(checkBuilt=TRUE,ask=FALSE)'
