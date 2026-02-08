#!/usr/bin/bash

# VirtualBox

sudo dnf update
sudo dnf install gcc kernel-devel kernel-headers dkms make bzip2 perl

sudo dnf install R
sudo dnf install R-devel
sudo dnf install bzip2-devel
sudo dnf install cairo-devel
sudo dnf install cmake
sudo dnf install freetype-devel
sudo dnf install fribidi-devel
sudo dnf install gcc-c++
sudo dnf install gcc-gfortran
sudo dnf install harfbuzz-devel
sudo dnf install java-21-openjdk-devel
sudo dnf install lapack-devel
sudo dnf install libcurl-devel
sudo yum install libdeflate-devel
sudo dnf install libjpeg-turbo-devel
sudo dnf install libjpeg-devel
sudo dnf install libpng-devel
sudo dnf install libtiff-devel
sudo dnf install libwebp-devel
sudo dnf install libxml2-devel
sudo dnf install libX11-devel
sudo dnf install libXt-devel
sudo dnf install libzstd-devel
sudo dnf install mesa-libGL-devel
sudo dnf install mysql-devel
sudo dnf install netcdf-devel
sudo dnf install pandoc
sudo dnf install pcre-devel
sudo dnf install pcre2-devel
sudo dnf install qpdf
sudo dnf install texlive
sudo dnf install fontconfig-devel
sudo dnf install readline-devel
sudo dnf install tcl tcl-devel tk tk-devel
sudo dnf install texlive-collection-latex
sudo dnf install texlive-collection-fontsextra
sudo dnf install texinfo
sudo dnf install texinfo-tex
sudo dnf install texlive-collection-fontsrecommended
sudo dnf install texlive-collection-latexrecommended
sudo dnf install tidy
sudo dnf install v8-devel
sudo dnf install xorg-x11-fonts-100dpi
sudo dnf install xorg-x11-fonts-75dpi
sudo dnf install xz-devel
# sudo dnf install java-1.8.0-openjdk-devel

# JAVA_HOME
export JAVAC=$(readlink -f $(which javac))
echo $JAVAC | sed 's|/bin/javac||'
export INCLUDE=$JAVA_HOME/include:$JAVA_HOME/include/linux:$INCLUDE
export PATH=$JAVA_HOME/bin:$PATH

cd ~
export R_LIBS=$HOME/R-devel/library
wget -qO- https://stat.ethz.ch/R/daily/R-devel.tar.gz | \
tar xvfz -
cd R-devel
./configure
make
ln -s ~/R-devel/bin/R ~/bin/R-devel
ln -s ~/R-devel/bin/Rscript ~/bin/Rscript-devel
Rscript-devel -e 'install.packages(c("shiny","V8"),repos="https://cran.r-project.org")'

# JAGS
export version=4.3.2
wget -qO- https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-${version}.tar.gz | \
tar xfz -
cd JAGS-${version}
./configure
make
sudo make install

# rjags
export PKG_CONFIG_PATH=/usr/local/lib/pkgconfig
export R_HOME=$HOME/R-devel
R-devel CMD INSTALL --configure-args='--enable-rpath' rjags

# library.zip

## Drop recommended packages
cat <<'EOL' | xargs -l -I {} zip -d library.zip {}*
base/
boot/
class/
cluster/
codetools/
compiler/
datasets/
foreign/
graphics/
grDevices/
grid/
KernSmooth/
lattice/
MASS/
Matrix/
methods/
mgcv/
nlme/
nnet/
parallel/
rpart/
spatial/
splines/
stats/
stats4/
survival/
tcltk/
tools/
translations/
utils/
EOL

## Add installed packages

cd R-devel/library
unzip ~/D/R/library

# ASAN

sudo dnf builddep R
sudo dnf install clang llvm gcc-c++ make texinfo
sudo dnf install llvm-symbolizer
svn checkout https://svn.r-project.org/R/trunk R-devel
cd R-devel
gfortran -print-file-name=libgfortran.so
gfortran -print-file-name=libquadmath.so
CC=clang \
CXX=clang++ \
FC=gfortran \
F77=gfortran \
FLIBS="-lgfortran -lquadmath" \
CFLAGS="-O1 -g -fsanitize=address,undefined" \
CXXFLAGS="-O1 -g -fsanitize=address,undefined" \
FFLAGS="-O1 -g" \
FCFLAGS="-O1 -g" \
LDFLAGS="-fsanitize=address,undefined" \
./configure \
  --prefix=$HOME/R-devel-asan \
  --enable-memory-profiling \
  --disable-java
export PATH="$HOME/R-devel-asan/bin:$PATH"
R --version

set -euo pipefail

# 1) Ensure ASAN runtime matches CRAN strictness
export ASAN_OPTIONS="detect_leaks=1:check_initialization_order=1:strict_init_order=1:halt_on_error=1"
export UBSAN_OPTIONS="halt_on_error=1"
export ASAN_SYMBOLIZER_PATH=/usr/bin/llvm-symbolizer

PKG_TARBALL=gap_1.10.tar.gz
echo "Running R CMD check --as-cran under ASAN..."
R CMD check "$PKG_TARBALL" --as-cran
