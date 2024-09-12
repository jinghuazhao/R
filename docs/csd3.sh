#!/usr/bin/bash

# 11/9/2024 onwards

module load pcre2/10.36/gcc/sya23vzi texlive/2015

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

export version=4.4.1
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
./configure --prefix=${prefix}/${dest} --with-pcre2 --enable-R-shlib
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
