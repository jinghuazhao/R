#!/usr/bin/env bash
set -euo pipefail

# CONFIGURATION
: "${R_SRC_DIR:=/home/jhz22/R-devel}"               # Path to R-devel source
: "${PREFIX:=/home/jhz22/R-asan}"                   # Installation directory
: "${JOBS:=$(nproc)}"                                # Parallel make jobs
: "${PKG_PATH:=/path/to/gap_1.14.tar.gz}"             # Path to the package
: "${R_LIBRARY_DIR1:=/home/jhz22/R-asan}"            # First package library path
: "${R_LIBRARY_DIR2:=/home/jhz22/R-devel/library}"   # Second package library path

# Sanitizer flags
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

# Set R library paths to include both directories
export R_LIBS_USER="${R_LIBRARY_DIR1}:${R_LIBRARY_DIR2}"

# Build R with ASAN
cd "${R_SRC_DIR}"
./configure --prefix="${PREFIX}" --enable-R-shlib --disable-openmp \
  CFLAGS="${CFLAGS}" CXXFLAGS="${CXXFLAGS}" FFLAGS="${FFLAGS}" LDFLAGS="${LDFLAGS}"
make -j"${JOBS}"
make install

# Ensure the new R binaries are used
export PATH="${PREFIX}/bin:${PATH}"

# Check the R installation
echo "=== Installed R version ==="
R --version

# Run install & check the package
echo "=== Running R CMD check --as-cran on gap_1.14.tar.gz ==="
R CMD check "${PKG_PATH}" --as-cran 2>&1 | tee check-as-cran.log

echo "=== DONE ==="
echo "Logs are in check-as-cran.log"
