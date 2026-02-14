#!/usr/bin/env bash

set -e

## ---- 1. Install ASan runtime (Fedora) ----
sudo dnf -y install libasan

## ---- 2. Set ASan compiler & linker flags ----
export CC="gcc -fsanitize=address,undefined -fno-omit-frame-pointer"
export CXX="g++ -fsanitize=address,undefined -fno-omit-frame-pointer"

export CFLAGS="-g -O2 -Wall -pedantic -fsanitize=address,undefined"
export CXXFLAGS="-g -O2 -Wall -pedantic -fsanitize=address,undefined"
export FCFLAGS="-g -O2"
export FFLAGS="-g -O2"

export MAIN_LDFLAGS="-fsanitize=address,undefined -pthread"

## Optional but often recommended
export ASAN_OPTIONS="detect_leaks=0:detect_odr_violation=0"

## ---- 3. Build R with ASan (if not already built) ----
## Adjust path to your R source directory
R_SRC=R-devel

if [ ! -d "$R_SRC" ]; then
  echo "R source directory $R_SRC not found."
  exit 1
fi

cd "$R_SRC"

./configure \
  --prefix=$HOME/R-asan \
  --enable-R-shlib \
  CC="$CC" CXX="$CXX" \
  CFLAGS="$CFLAGS" CXXFLAGS="$CXXFLAGS" \
  MAIN_LDFLAGS="$MAIN_LDFLAGS"

make -j$(nproc)
make install

cd ..

## ---- 4. Run package check with ASan R ----
export PATH=$HOME/R-asan/bin:$PATH

R CMD check --as-cran gap_1.14.tar.gz
