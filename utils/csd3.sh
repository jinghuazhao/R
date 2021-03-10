#!/usr/bin/bash

module load gcc/6
module load pcre/8.38
module load texlive
export version=4.0.4

IFS=\. read -a fields <<<${version}
export major=${fields[0]}
export minor1=${fields[1]}
export minor2=${fields[2]}
echo ${major}.${minor1}.${minor2}

wget https://cran.r-project.org/src/base/R-${main}/R-${version}.tar.gz
tar xvfz R-${version}.tar.gz
cd ${version}
export prefix=/rds-d4/user/$USER/hpc-work
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${prefix}/include LDFLAGS=-L${prefix}/lib
make
make install
cd $HOME/bin
ln -sf  /rds-d4/user/$USER/hpc-work/R-${version}/bin/R
