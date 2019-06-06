# utilities

The standard CRAN submission check is via [check.sh](check.sh).

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

A more sophisticated way is as in [install.sh](install.sh).

This should be possible with configure as well, i.e., [configure.sh](configure.sh), but I reckon this is hardly necessary.
