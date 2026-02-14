#!/usr/bin/bash

export PKG_TARBALL=gap_1.14.tar.gz

function asan_gcc()
# gcc
{
export R_LIBS=$HOME/R-devel-gcc/library:$HOME/R-devel/library
export PATH=$HOME/R-devel-gcc/bin:$PATH
R CMD check $PKG_TARBALL --as-cran
}

function asan_llvm()
{
# clang
export R_LIBS=$HOME/R-devel-clang/library:$HOME/R-devel/library
export PATH=$HOME/R-devel-clang/bin:$PATH
ASAN_OPTIONS=detect_leaks=0 \
UBSAN_OPTIONS=print_stacktrace=1 \
R CMD check $PKG_TARBALL \
  --no-manual
}

function asan()
# asan
{
set -euo pipefail
export ASAN_OPTIONS="detect_leaks=1:check_initialization_order=1:strict_init_order=1:halt_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1"
export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer
R CMD check "$PKG_TARBALL" --as-cran
}

asan

function compileR()
{
cd R-devel

## GCC
export CC=gcc
export CXX=g++
export FC=gfortran
export F77=gfortran

export CFLAGS="-O2 -g"
export CXXFLAGS="-O2 -g"
export FFLAGS="-O2 -g"
export FCFLAGS="-O2 -g"

unset LDFLAGS
unset FLIBS

./configure \
  --prefix=$HOME/R-devel-gcc \
  --enable-memory-profiling \
  --disable-java
make -j1
make install

## clang
export CC=clang
export CXX=clang++
export FC=flang
export F77=flang

export CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"
export CXXFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"
export FFLAGS="-O1 -g"
export FCFLAGS="-O1 -g"
export LDFLAGS="-fsanitize=address,undefined"

unset FLIBS

./configure \
  --prefix=$HOME/R-devel-llvm \
  --enable-memory-profiling \
  --disable-java
make -j1
make install

## ASAN (CRAN-style clang + gfortran)
gfortran -print-file-name=libgfortran.so
gfortran -print-file-name=libquadmath.so
export CC=clang
export CXX=clang++
export FC=gfortran
export F77=gfortran
export FLIBS="-lgfortran -lquadmath"
export CFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"
export CXXFLAGS="-O1 -g -fsanitize=address,undefined -fno-omit-frame-pointer"
export FFLAGS="-O1 -g"
export FCFLAGS="-O1 -g"
export LDFLAGS="-fsanitize=address,undefined"
./configure \
  --prefix=$HOME/R-devel-asan \
  --enable-memory-profiling \
  --disable-java
make -j1
make install
export LD_LIBRARY_PATH=$(dirname $(gfortran -print-file-name=libgfortran.so)):$LD_LIBRARY_PATH
export R_LIBS=$HOME/R-devel-asan/library:$HOME/R-devel/library
export PATH="$HOME/R-devel-asan/bin:$PATH"
export R_ENVIRON_USER=
export R_PROFILE_USER=
R --version
ASAN_OPTIONS=detect_leaks=0:alloc_dealloc_mismatch=0 \
UBSAN_OPTIONS=print_stacktrace=1 \
R CMD check $PKG_TARBALL \
  --as-cran \
  --no-manual
}
