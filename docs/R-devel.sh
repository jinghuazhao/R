#!/usr/bin/bash

sudo dnf install R
sudo dnf install R-devel
sudo dnf install gcc-c++
sudo dnf install gcc-gfortran
sudo dnf install pcre-devel
sudo dnf install pcre2-devel
sudo dnf install java-1.8.0-openjdk-devel
sudo dnf install readline-devel
sudo dnf install libcurl-devel
sudo dnf install libjpeg-turbo-devel
sudo dnf install libX11-devel
sudo dnf install libXt-devel
sudo dnf install bzip2-devel
sudo dnf install xz-devel
sudo dnf install pandoc
sudo dnf install qpdf
sudo dnf install texlive-collection-latex
sudo dnf install texlive-collection-fontsextra
sudo dnf install texinfo
sudo dnf install texinfo-tex
sudo dnf install texlive-collection-fontsrecommended
sudo dnf install texlive-collection-latexrecommended
sudo dnf install tidy
sudo dnf install lapack-devel
sudo dnf install v8-devel
sudo dnf install xorg-x11-fonts-100dpi
sudo dnf install xorg-x11-fonts-75dpi

cd ~
export R_LIBS=$HOME/R-devel/library
wget -qO- https://stat.ethz.ch/R/daily/R-devel.tar.gz | \
tar xvfz -
cd R-devel
./configure
make

# symbolic links

ln -s ~/R-devel/bin/R ~/bin/R-devel
ln -s ~/R-devel/bin/Rscript ~/bin/Rscript-devel
Rscript-devel -e 'install.packages(c("shiny","V8"),repos="https://cran.r-project.org")'

# JAGS

export version=4.3.2
wget -qO- https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-${version}.tar.gz | \
tar xfz -
cd JAGS-${version}
./configure --without-recommended-packages
make
sudo make install

# Added packages from transferred R-devel/library

cat <<'EOL' | xargs -l -I {} zip -d library.zip R-devel/library/{}
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
valr/
EOL
