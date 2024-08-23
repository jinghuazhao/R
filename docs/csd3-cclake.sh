#!/usr/bin/bash

module load gcc/6
module load curl-7.63.0-gcc-5.4.0-4uswlql geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc
module load mono/5.0.1.1 netcdf/4.4.1
module load pcre2-10.20-gcc-5.4.0-tcuhtrb protobuf-3.4.0-gcc-5.4.0-zkpendv texlive
module load ceuadmin/brotli/1.0.9 ceuadmin/icu/70.1 ceuadmin/jq/1.6
module load ceuadmin/glpk/4.57 ceuadmin/NLopt/2.7.1
module load ceuadmin/readline/8.0 ceuadmin/rtmpdump/2.3

export version=4.4.1
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
bin/Rscript -e 'update.packages(checkBuilt=TRUE,ask=FALSE)'

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
