#!/usr/bin/bash

CFLAGS="-g -O2 -Wall -pedantic -mtune=native -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong -fstack-clash-protection -fcf-protection -Werror=implicit-function-declaration -Wstrict-prototypes" \
FFLAGS="-g -O2 -mtune=native -Wall -pedantic" \
CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-parentheses -Werror=format-security -Wp,-D_FORTIFY_SOURCE=2 -fexceptions -fstack-protector-strong -fstack-clash-protection -fcf-protection" \
JAVA_HOME=/usr/lib/jvm/java-11 \
./configure
make
R-devel CMD check --as-cran gap_1.3.tar.gz

# https://www.stats.ox.ac.uk/pub/bdr/Rconfig/r-devel-linux-x86_64-fedora-gcc
