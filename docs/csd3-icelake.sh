#!/usr/bin/bash

module load gcc/8 geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
module load curl/7.83.0/gcc/ozlrq5hx
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc readline/8.1/gcc/bgw44yb2
module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1 ceuadmin/nettle/2.7.1

export version=4.3.3
IFS=\. read major minor1 minor2 <<<${version}
export prefix=/rds/project/jmmh2/rds-jmmh2-public_databases/software/R-${version}-icelake
mkdir -p ${prefix} && cd ${prefix}
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz - -C ${prefix}
mv */* .
rmdir R-${version}
export gcc8=/usr/local/software/master/gcc/8
export intl=/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/gettext-0.19.8.1-5iqkv5zrractwd57vydu5czosgtrlwj2/
export HPC_WORK=/rds/user/jhz22/hpc-work
export include=${gcc8}/include:${intl}/include:${HPC_WORK}/include
export ldflags=${gcc8}/lib64:${gcc8}/lib:${intl}/lib:${HPC_WORK}/lib64:${HPC_WORK}/lib
./configure --prefix=${prefix}/R-${version}-icelake \
            --with-pcre2 \
            --enable-R-shlib CPPFLAGS=-I${include} LDFLAGS=-L${ldflags} LIBS=-ltinfo LIBS=-lintl
make
