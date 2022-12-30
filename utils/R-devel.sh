#!/usr/bin/bash

sudo dnf install R
sudo dnf install R-devel
sudo dnf install gcc-c++
sudo dnf install gcc-gfortran
sudo dnf install pcre-devel
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
sudo dnf install texinfo-tex
sudo dnf install texlive-collection-fontsrecommended
sudo dnf install texlive-collection-latexrecommended
sudo dnf install tidy
sudo dnf install lapack-devel
sudo dnf install v8-devel
sudo dnf install xorg-x11-fonts-100dpi
sudo dnf install xorg-x11-fonts-75dpi
export R_LIBS=$HOME/R-devel/library
wget -qO- https://stat.ethz.ch/R/daily/R-devel.tar.gz | \
tar xvfz -
cd R-devel
./configure
make
wget -qO- https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Source/JAGS-4.3.1.tar.gz | \
tar xfz -
./configure
make
sudo make install
