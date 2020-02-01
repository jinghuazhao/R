## utilities

* [check.sh](check.sh) is the standard CRAN submission check, e.g.,
```bash
check.sh gap
```
* [install.sh](install.sh) is a sophistic way of installation.
One could implement a more restrictive gcc9 specifications, as described in [https://www.stats.ox.ac.uk/pub/bdr/gcc9/README.txt](https://www.stats.ox.ac.uk/pub/bdr/gcc9/README.txt), as follows,
```bash
export CC=/usr/bin/gcc
export CXX=/usr/g++
export FC=/usr/bin/gfortran
export CFLAGS="-g -O2 -Wall -pedantic -mtune=native"
export FFLAGS="-g -O2 -mtune=native -Wall -pedantic"
export CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses"
export LDFLAGS="-L/usr/lib64 -L/usr/lib64"
```
followed by `R-devel CMD INSTALL gap_1.2.2.tar.gz`, say.

The CPPFLAGS="-D_FORTIFY_SOURCE=2" produces

> Warning: ignoring return value of ‘fgets’, declared with attribute warn_unused_result [-Wunused-result]

* [configure.sh](configure.sh) is notable for compiling new release of R.

* [register.R](register.R) can be used as 
```bash
[Rscript] register.R gap
```
to produce `package_native_routine_registration_skeleton.c`.

* [st.sh](st.sh) is a batch file for GitHub.

* [gap-suggests.R](gap-suggests.R) is useful to install susggested packages for gap.
