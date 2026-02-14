#!/usr/bin/env bash

: "${R_SRC_DIR:=/home/jhz22/R-devel}"                # Path to R-devel source
: "${PREFIX:=/home/jhz22/R-asan}"                    # Installation directory
: "${JOBS:=$(nproc)}"                                # Parallel make jobs
: "${PKG_PATH:=/path/to/gap_1.14.tar.gz}"            # Path to the package
: "${R_LIBRARY_DIR1:=/home/jhz22/R-asan}"            # First package library path
: "${R_LIBRARY_DIR2:=/home/jhz22/R-devel/library}"   # Second package library path

ASAN_FLAGS="-fsanitize=address,undefined,bounds-strict"
COMMON_CFLAGS="-g -O1 -Wall -Wextra"

export CC="clang"
export CXX="clang++"
export FC="flang"
export AR="llvm-ar"
export RANLIB="llvm-ranlib"
export NM="llvm-nm"
export LD="clang++"

export CFLAGS="${COMMON_CFLAGS} ${ASAN_FLAGS}"
export CXXFLAGS="${COMMON_CFLAGS} ${ASAN_FLAGS} -frtti"
export FFLAGS="-g -O2 -Wall"
export LDFLAGS="${ASAN_FLAGS}"

export R_LIBS_USER="${R_LIBRARY_DIR1}:${R_LIBRARY_DIR2}"

cd "${R_SRC_DIR}"
./configure --prefix="${PREFIX}" --enable-R-shlib --disable-openmp \
  CFLAGS="${CFLAGS}" CXXFLAGS="${CXXFLAGS}" FFLAGS="${FFLAGS}" LDFLAGS="${LDFLAGS}"
make -j"${JOBS}"
make install

export PATH="${PREFIX}/bin:${PATH}"
R --version
R CMD check "${PKG_PATH}" --as-cran 2>&1 | tee check-as-cran.log
