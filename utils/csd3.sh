#!/usr/bin/bash

module load gcc/6
module load pcre/8.38
module load texlive

export version=4.0.5
export major=$(cut -d. -f1 <<<${version})
export minor1=$(cut -d. -f2 <<<${version})
export minor2=$(cut -d. -f3 <<<${version})
echo ${major}.${minor1}.${minor2}

export prefix=/rds-d4/user/$USER/hpc-work
cd ${prefix}
wget https://cran.r-project.org/src/base/R-${major}/R-${version}.tar.gz
tar xvfz R-${version}.tar.gz
cd R-${version}
./configure --prefix=${prefix} \
            --with-pcre1 \
            --enable-R-shlib CPPFLAGS=-I${prefix}/include LDFLAGS=-L${prefix}/lib
make
make install
cd $HOME/bin
ln -sf  ${prefix}/R-${version}/bin/R

# ---

function read_parse_version()
{
  IFS=\. read -a fields <<<${version}
  export major=${fields[0]}
  export minor1=${fields[1]}
  export minor2=${fields[2]}
}

