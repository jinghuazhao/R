#!/usr/bin/bash

module load gcc/9 geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc
module load ceuadmin/brotli/1.0.9 ceuadmin/glpk/4.57 ceuadmin/icu/70.1 ceuadmin/jq/1.6 protobuf-3.4.0-gcc-5.4.0-zkpendv
module load ceuadmin/NLopt/2.7.1

export version=4.4.0
IFS=\. read major minor1 minor2 <<<${version}
export rds=/rds/project/jmmh2/rds-jmmh2-public_databases/software
export prefix=$CEUADMIN/R
export dest=${version}
export R_LIBS=${rds}/R
cd ${prefix}
mkdir ${dest}
umask 022
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xfz - --no-same-owner -C ${dest} --strip-components=1
cd ${dest}
./configure --prefix=${prefix}/${dest} --with-pcre2 --enable-R-shlib
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
