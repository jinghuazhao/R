#!/usr/bin/env bash
set -euo pipefail
cd ~/R

export _R_CHECK_CRAN_INCOMING_=TRUE
export _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE
export _R_CHECK_FORCE_SUGGESTS_=FALSE
export _R_CHECK_LIMIT_CORES_=TRUE
export _R_CHECK_BUILD_VIGNETTES_=FALSE
export MAKEFLAGS="-j$(nproc)"

echo "Using R:"
Rscript-devel -e 'cat(R.version.string,"\n")'

# 1) Update Rd files
Rscript-devel -e 'devtools::document("gap")'

# 2) Build package (THIS builds and indexes vignettes)
R-devel CMD build --compact-vignettes=both gap

# 3) Install & check exactly like CRAN
TARBALL=$(ls gap_*.tar.gz | tail -1)
R-devel CMD INSTALL "$TARBALL"
R-devel CMD check --as-cran --run-donttest "$TARBALL"

rm -f gap/src/*.so gap/src/*.o || true
