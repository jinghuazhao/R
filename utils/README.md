## utilities

This directory contains various utility scripts.

* [check.sh](check.sh) is the standard CRAN submission check, e.g.,
```bash
check.sh gap
```
* [csd3.sh](csd3.sh) contains steps to install R under CSD3. Note that for `check.sh` work properly, it is necessary to load the modules in it.

* [install.sh](install.sh) is a sophisticated way of installation. One could implement a more restrictive gcc9 specifications, as described in [https://www.stats.ox.ac.uk/pub/bdr/gcc9/README.txt](https://www.stats.ox.ac.uk/pub/bdr/gcc9/README.txt), as follows,
```bash
export CC=/usr/bin/gcc
export CXX=/usr/g++
export FC=/usr/bin/gfortran
export CFLAGS="-g -O2 -Wall -pedantic -mtune=native"
export FFLAGS="-g -O2 -mtune=native -Wall -pedantic"
export CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses"
export LDFLAGS="-L/usr/lib64 -L/usr/lib64"
```
followed by `R-devel CMD INSTALL gap_1.2.3.tar.gz`, say. They could also be part of $HOME/.R/Makevars, omitting `export`.

The CPPFLAGS="-D_FORTIFY_SOURCE=2" produces

> Warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result [-Wunused-result]

It also shows how to install `gap` from GitHub via `Rscript`.

* [register.R](register.R) can be used as 
```bash
[Rscript] register.R gap
```
to produce `package_native_routine_registration_skeleton.c`.

It is more apporopriate for `gfortran` 9.2 and later to extract C prototypes for Fortran subroutines with a special flag:
```bash
gfortran -c -fc-prototypes-external lmm.f
```
as in package 'lmm'.

* [st.sh](st.sh) is a batch file for GitHub.

* [suggests.R](suggests.R) is useful to install susggested packages.
