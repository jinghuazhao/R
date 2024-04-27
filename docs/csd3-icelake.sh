#!/usr/bin/bash

module load curl/7.79.0/gcc/75dxv7ac pcre2/10.36/gcc/sya23vzi texlive
module load ceuadmin/json-c/0.17-20230812-icelake ceuadmin/openssl/3.2.1-icelake
module load ceuadmin/nettle/2.7.1

export version=4.4.0
IFS=\. read major minor1 minor2 <<<${version}
export prefix=/rds/project/jmmh2/rds-jmmh2-public_databases/software
export R_LIBS=${prefix}/R-icelake
cd ${prefix}
mkdir R-${version}-icelake
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz - -C R-${version}-icelake --strip-components=1
cd R-${version}-icelake
./configure --prefix=${prefix}/R-${version}-icelake --with-pcre2 --enable-R-shlib
make
make install

function legacy()
{
  module load gcc/8 geos-3.6.2-gcc-5.4.0-vejexvy libiconv/1.16/gcc/4miyzf3w pcre2/10.36/gcc/sya23vzi texlive
  module load curl/7.83.0/gcc/ozlrq5hx gettext/0.21/gcc/qnrcglqo
  module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc libpng/1.6.37/intel/jfrl6z6c readline/8.1/gcc/bgw44yb2
  module load ceuadmin/json-c/0.17-20230812-icelake
  module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1 ceuadmin/nettle/2.7.1 ceuadmin/openssl/3.2.1-icelake
}
