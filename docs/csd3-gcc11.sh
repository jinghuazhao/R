#!/usr/bin/bash

module load gcc/11 geos-3.6.2-gcc-5.4.0-vejexvy gettext-0.19.8.1-gcc-5.4.0-5iqkv5z pcre2-10.20-gcc-5.4.0-tcuhtrb texlive
module load image-magick-7.0.5-9-gcc-5.4.0-d4lemcc
module load ceuadmin/glpk/4.57 ceuadmin/icu/70.1

export version=4.3.3
export prefix=/rds/project/jmmh2/rds-jmmh2-public_databases/software/R-${version}-gcc11
cd ${prefix}
./configure --prefix=${prefix} \
            --with-pcre2 \
            --enable-R-shlib
make
