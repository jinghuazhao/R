#!/usr/bin/bash

module load curl/7.83.0/gcc/ozlrq5hx
module load gcc/11.3.0/gcc/4zpip55j
module load geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc libpng/1.6.37/gcc/bkdpz5q4 protobuf-3.4.0-gcc-5.4.0-zkpendv
module load ceuadmin/brotli/1.1.0 ceuadmin/glpk/4.57 ceuadmin/icu/70.1 ceuadmin/jq/1.6 ceuadmin/NLopt/2.7.1
module load curl/7.79.0/gcc/75dxv7ac ceuadmin/libiconv/1.17

export version=4.5.1
IFS=\. read major minor1 minor2 <<<${version}
export rds=/rds/project/jmmh2/rds-jmmh2-public_databases/software
export prefix=${CEUADMIN}/R
export dest=${version}-gcc11
export R_LIBS=${rds}/R-gcc11
cd ${prefix}
mkdir ${dest}
umask 022
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xfz - --no-same-owner -C ${dest} --strip-components=1
cd ${dest}
./configure --prefix=${prefix}/${dest} --with-pcre2 --enable-R-shlib \
            CPPFLAGS=-I/usr/local/Cluster-Apps/ceuadmin/libiconv/1.17/include \
            LDFLAGS='-L/usr/local/Cluster-Apps/ceuadmin/libiconv/1.17/lib -liconv'
make
