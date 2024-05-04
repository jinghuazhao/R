#!/usr/bin/bash

module load curl/7.79.0/gcc/75dxv7ac gettext/0.21/gcc/qnrcglqo libiconv/1.16/intel/64iicvbf
module load libpng/1.6.37/intel/jfrl6z6c pcre2/10.36/gcc/sya23vzi readline/8.1/gcc/bumlt4j6
module load ceuadmin/json-c/0.17-20230812-icelake
module load ceuadmin/nettle/2.7.1 texlive

export version=4.4.0
IFS=\. read major minor1 minor2 <<<${version}
export rds=/rds/project/jmmh2/rds-jmmh2-public_databases/software
export prefix=$CEUADMIN/R
export dest=${version}-icelake
export R_LIBS=${rds}/R-icelake
cd ${prefix}
mkdir ${dest}
umask 022
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar --no-same-owner xfz - -C ${dest} --strip-components=1
cd ${dest}
./configure --prefix=${prefix}/${dest} --with-pcre2 --enable-R-shlib
make

function legacy()
{
  module load gcc/8 geos-3.6.2-gcc-5.4.0-vejexvy libiconv/1.16/gcc/4miyzf3w pcre2/10.36/gcc/sya23vzi texlive
  module load curl/7.83.0/gcc/ozlrq5hx gettext/0.21/gcc/qnrcglqo
  module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc libpng/1.6.37/intel/jfrl6z6c readline/8.1/gcc/bgw44yb2
  module load ceuadmin/json-c/0.17-20230812-icelake
  module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1 ceuadmin/nettle/2.7.1 ceuadmin/openssl/3.2.1-icelake
}
