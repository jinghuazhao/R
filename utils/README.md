## utilities

This directory contains various utility scripts.

* [check.sh](check.sh) is the standard CRAN submission check, e.g.,
```bash
check.sh gap
```
* [csd3.sh](csd3.sh) contains steps to install R under CSD3. Note that for `check.sh` to work properly, it is necessary to load the modules in it.

* [install.sh](install.sh) is a sophisticated way of installation. One could implement a more restrictive gcc9 specifications, as orignally described in [https://www.stats.ox.ac.uk/pub/bdr/gcc9/README.txt](https://www.stats.ox.ac.uk/pub/bdr/gcc9/README.txt), as follows,
```bash
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export FC=/usr/bin/gfortran
export CFLAGS="-g -O2 -Wall -pedantic -mtune=native"
export FFLAGS="-g -O2 -mtune=native -Wall -pedantic"
export CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses -Warray-parameter"
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

pQTLtools is somewhat poorly mirrored here, [https://rdrr.io/github/jinghuazhao/pQTLtools/](https://rdrr.io/github/jinghuazhao/pQTLtools/).

## clang / gcc errror messages

The settings below are derived from the following links,

* clang, <https://www.stats.ox.ac.uk/pub/bdr/Rconfig/r-devel-linux-x86_64-fedora-clang>
* gcc, <https://www.stats.ox.ac.uk/pub/bdr/Rconfig/r-devel-linux-x86_64-fedora-gcc>

Both ininove Fedora 36, which has gcc 12.0.1 and clang 14.0.5 although the error mesages were indicated as for clang 15 on CRAN.

Fedora 36 setup gets easier to start with 'sudo dnf install R-devel' followed by adding `cmake`, `pandoc`, `ImageMagick` as well as `cairo-devel`, `libcurl-devel`, `libjpeg-turbo-devel`, `readline-devel`, `v8-devel`, `xorg-x11-fonts*`.

Note that some errors can only be seen through R CMD INSTALL.

### .R/Makevars

```bash
export _R_CHECK_INSTALL_DEPENDS_ true
## the next is the default, but --as-cran has true
export _R_CHECK_SUGGESTS_ONLY_ false
export _R_CHECK_NO_RECOMMENDED_ true
export _R_CHECK_DOC_SIZES2_ true
export _R_CHECK_DEPRECATED_DEFUNCT_ true
export _R_CHECK_SCREEN_DEVICE_ warn
export _R_CHECK_REPLACING_IMPORTS_ true
export _R_CHECK_TOPLEVEL_FILES_ true
export _R_CHECK_DOT_FIRSTLIB_ true
export _R_CHECK_RD_LINE_WIDTHS_ true
export _R_CHECK_S3_METHODS_NOT_REGISTERED_ true
export _R_CHECK_OVERWRITE_REGISTERED_S3_METHODS_ true
export _R_CHECK_CODE_USAGE_WITH_ONLY_BASE_ATTACHED_ TRUE
export _R_CHECK_NATIVE_ROUTINE_REGISTRATION_ true
export _R_CHECK_FF_CALLS_ registration
export _R_CHECK_PRAGMAS_ true
export _R_CHECK_COMPILATION_FLAGS_ true
export _R_CHECK_R_DEPENDS_ true
export _R_CHECK_PACKAGES_USED_IN_TESTS_USE_SUBDIRS_ true
export _R_CHECK_PKG_SIZES_ false
export _R_CHECK_SHLIB_OPENMP_FLAGS_ true

export _R_CHECK_LIMIT_CORES_ true
export _R_CHECK_LENGTH_1_CONDITION_ package:_R_CHECK_PACKAGE_NAME_,verbose
#export _R_CHECK_LENGTH_1_LOGIC2_ "package:_R_CHECK_PACKAGE_NAME_,verbose"
export _R_S3_METHOD_LOOKUP_BASEENV_AFTER_GLOBALENV_ true
export _R_CHECK_COMPILATION_FLAGS_KNOWN_ "-Wno-deprecated-declarations -Wno-ignored-attributes -Wno-parentheses-Werror=format-security -Wp,-D_FORTIFY_SOURCE=2i -Werror=implicit-function-declaration"
export _R_CHECK_AUTOCONF_ true
export _R_CHECK_THINGS_IN_CHECK_DIR_ true
export _R_CHECK_THINGS_IN_TEMP_DIR_ true
export _R_CHECK_THINGS_IN_TEMP_DIR_EXCLUDE_ "^ompi"
export _R_CHECK_BASHISMS_ true
export _R_CHECK_DEPENDS_ONLY_DATA_ true
export _R_CHECK_MATRIX_DATA_ TRUE
export _R_CHECK_RD_VALIDATE_RD2HTML_ true
export _R_CHECK_RD_MATH_RENDERING_ true
```

### R-devel

This is how R-devel is compiled,

```bash
cd ${HOME}
wget -qO- https://cran.r-project.org/src/base-prerelease/R-devel.tar.gz | \
tar xfz -
cd R-devel

CFLAGS="-g -O2 -Wall -pedantic -mtune=native -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong -fstack-clash-protection -fcf-protection -Werror=implicit-function-declaration -Wstrict-prototypes" \
FFLAGS="-g -O2 -mtune=native -Wall -pedantic" \
CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-parentheses -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong -fstack-clash-protection -fcf-protection" \
JAVA_HOME=/usr/lib/jvm/java-11 \
./configure
make
ln -sf R-devel/bin/R ${HOME}/bin/R-devel
```

A version of these is noted in [R-devel.sh](R-devel.sh) for Fedora 37.

### R CMD check

This is standard,

```bash
R-devel CMD check --as-cran gap_1.3.tar.gz
```

For a while under Fedora 36, there has been the following error message,

```
* checking re-building of vignette outputs ... [34s/45s] ERROR
Error(s) in re-building vignettes:
  ...
--- re-building ‘gap.Rmd’ using rmarkdown
Read 16 records
Warning in grid.Call.graphics(C_segments, x$x0, x$y0, x$x1, x$y1, x$arrow) :
  semi-transparency is not supported on this device: reported only once per page
Failed with error:  'error reading from connection'
Quitting from lines 668-670 (gap.Rmd) 
Error: processing vignette 'gap.Rmd' failed with diagnostics:
error reading from connection
--- failed re-building ‘gap.Rmd’

--- re-building ‘shinygap.Rmd’ using rmarkdown
--- finished re-building ‘shinygap.Rmd’

--- re-building ‘jss.Rnw’ using Sweave
--- finished re-building ‘jss.Rnw’

SUMMARY: processing the following file failed:
  ‘gap.Rmd’

Error: Vignette re-building failed.
Execution halted
```

It turns out package `meta` is missing from the package list, which removes the error after installation.

It is known that some preparations are needed to compile PDF files, e.g., 

```
* checking PDF version of manual ... WARNING
LaTeX errors when creating PDF version.
This typically indicates Rd problems.
LaTeX errors found:
! LaTeX Error: File `inconsolata.sty' not found.

Type X to quit or <RETURN> to proceed,
or enter new name. (Default extension: sty)

! Emergency stop.
<read *>

l.303 ^^M

!  ==> Fatal error occurred, no output PDF file produced!
* checking PDF version of manual without index ... OK
* DONE

Status: 1 WARNING, 1 NOTE
See
  ‘/home/jhz22/R/gap.examples.Rcheck/00check.log’
for details.
```

On Cambridge University's CSD3 system, one only needs to use `module load texlive` first.

## html

This is standard, e.g.,

```bash
pandoc README.md --citeproc --mathjax -s --self-contained -o index.html
```

## R version

R version x.x.x. can be parsed as follows,

```bash
export version=4.2.2
IFS=\. read major minor1 minor2 <<<${version}
echo ${major}.${minor1}.${minor2}
```
