#!/usr/bin/bash

R-devel CMD check --as-cran $1

export CC=/usr/bin/gcc
export CXX=/usr/g++
export FC=/usr/bin/gfortran
export CFLAGS="-g -O2 -Wall -pedantic -mtune=native"
export FFLAGS="-g -O2 -mtune=native -Wall -pedantic"
export CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses"
export LDFLAGS="-L/usr/lib64 -L/usr/lib64"

R-devel CMD INSTALL $1
