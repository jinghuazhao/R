#!/usr/bin/bash

# 31/10/2024 onwards
module load curl/7.83.0/gcc/ozlrq5hx
module load hdf5/1.12.1 icu4c/67.1/gcc/maavowaj libpng/1.6.37/intel/jfrl6z6c
module load mono/5.0.1.1 netcdf/4.4.1
module load pcre2/10.36/gcc/sya23vzi texlive/2015
module load jags-4.3.0-gcc-5.4.0-4z5shby
module load ceuadmin/libsodium ceuadmin/rust ceuadmin/libiconv/1.17 ceuadmin/NLopt/2.7.1
module load rstudio/2024.04.2+764

function my_load()
{
  module load curl/7.79.0/gcc/75dxv7ac gettext/0.21/gcc/qnrcglqo libiconv/1.16/intel/64iicvbf
  module load libpng/1.6.37/intel/jfrl6z6c pcre2/10.36/gcc/sya23vzi readline/8.1/gcc/bumlt4j6
  module load texlive
  module load libdeflate/1.10/gcc/6ij3yqv2
  module load ceuadmin/json-c/0.17-20230812-icelake ceuadmin/krb5/1.21.2-icelake
  module load ceuadmin/nettle/3.9-icelake ceuadmin/qpdf/11.9.1

  module unload gcc/6
}

export version=4.5.0
IFS=\. read major minor1 minor2 <<<${version}
export rds=/rds/project/rds-4o5vpvAowP0/software
export prefix=$CEUADMIN/R
export dest=${version}-icelake
umask 022
cd ${prefix}
mkdir ${dest}
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xfz - --no-same-owner -C ${dest} --strip-components=1
cd ${dest}
./configure --prefix=${prefix}/${dest} --with-pcre2 --enable-R-shlib \
            CPPFLAGS=-I/usr/local/Cluster-Apps/ceuadmin/libiconv/1.17/include \
            LDFLAGS='-L/usr/local/Cluster-Apps/ceuadmin/libiconv/1.17/lib -liconv'
make
export R_LIBS=${rds}/R

function legacy()
{
  module load gcc/8 geos-3.6.2-gcc-5.4.0-vejexvy libiconv/1.16/gcc/4miyzf3w pcre2/10.36/gcc/sya23vzi texlive
  module load curl/7.83.0/gcc/ozlrq5hx gettext/0.21/gcc/qnrcglqo
  module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc libpng/1.6.37/intel/jfrl6z6c readline/8.1/gcc/bgw44yb2
  module load ceuadmin/json-c/0.17-20230812-icelake
  module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1 ceuadmin/nettle/2.7.1 ceuadmin/openssl/3.2.1-icelake
}
