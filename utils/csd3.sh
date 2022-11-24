#!/usr/bin/bash

module load gcc/6 texlive

export version=4.2.2
export major=$(cut -d. -f1 <<<${version})
export minor1=$(cut -d. -f2 <<<${version})
export minor2=$(cut -d. -f3 <<<${version})
echo ${major}.${minor1}.${minor2}

export prefix=/rds-d4/user/$USER/hpc-work
cd ${prefix}
wget -qO- https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz | \
tar xvfz -
cd R-${version}
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${prefix}/include LDFLAGS=-L${prefix}/lib
make
make install
cd $HOME/bin
ln -sf  ${prefix}/R-${version}/bin/R
Rscript -e 'update.packages(checkBuilt=TRUE,ask=FALSE)'

# --- more recent pcre has been installed independently

module load pcre/8.38

function read_parse_version()
{
  IFS=\. read -a fields <<<${version}
  export major=${fields[0]}
  export minor1=${fields[1]}
  export minor2=${fields[2]}
}

