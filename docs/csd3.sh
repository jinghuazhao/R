#!/usr/bin/bash

module load gcc/6 geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc
module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1

export prefix=/rds/project/jmmh2/rds-jmmh2-public_databases/software
export R_LIBS=${prefix}/R
cd ${prefix}
export version=4.4.0
IFS=\. read major minor1 minor2 <<<${version}
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz -
cd R-${version}
./configure --prefix=${prefix}/R-${version} --with-pcre2 --enable-R-shlib
make
make install
Rscript -e 'update.packages(checkBuilt=TRUE,ask=FALSE)'

function pre_4.3.3()
{
  export gcc6=/usr/local/software/master/gcc/6
  export intl=/usr/local/software/spack/spack-0.11.2/opt/spack/linux-rhel7-x86_64/gcc-5.4.0/gettext-0.19.8.1-5iqkv5zrractwd57vydu5czosgtrlwj2/
  export HPC_WORK=/rds/user/jhz22/hpc-work
  export include=${gcc6}/include:${intl}/include:${HPC_WORK}/include
  export ldflags=${gcc6}/lib64:${gcc6}/lib:${intl}/lib:${HPC_WORK}/lib64:${HPC_WORK}/lib
  cd R-${version}
  ./configure --prefix=${prefix} \
              --with-pcre2 \
              --enable-R-shlib CPPFLAGS=-I${include} LDFLAGS=-L${ldflags} LIBS=-ltinfo LIBS=-lintl
}
