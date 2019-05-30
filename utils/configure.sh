#!/usr/bin/bash

./configure --enable-R-shlib \
 CC=/usr/bin/gcc \
 CXX=/usr/g++ \
 FC=/usr/bin/gfortran \
 CFLAGS="-g -O2 -Wall -pedantic -mtune=native" \
 FFLAGS="-g -O2 -mtune=native -Wall -pedantic" \
 CXXFLAGS="-g -O2 -Wall -pedantic -mtune=native -Wno-ignored-attributes -Wno-deprecated-declarations -Wno-parentheses" \
 LDFLAGS="-L/usr/lib64 -L/usr/lib64"
