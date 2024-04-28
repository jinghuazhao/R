#!/usr/bin/bash

module load gcc/11 geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc
module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1

export version=4.4.0
IFS=\. read major minor1 minor2 <<<${version}
export prefix=/rds/project/jmmh2/rds-jmmh2-public_databases/software
export R_LIBS=${prefix}/R-gcc11
cd ${prefix}
mkdir R-${version}-gcc11
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz - -C R-${version}-gcc11 --strip-components=1
cd R-${version}-gcc11
export prefix=/usr/local/Cluster-Apps/ceuadmin/R
./configure --prefix=${prefix}/${version}-gcc11 --with-pcre2 --enable-R-shlib
make
make install
