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

Rscript-devel -e '
pkg <- "gap"

cat("Using R:", R.version.string, "\n")

# 1) generate Rd
devtools::document(pkg)

# 2) build tarball (modern replacement for R CMD build)
tarball <- pkgbuild::build(pkg, clean = TRUE)

# 3) install from source
install.packages(tarball, repos = NULL, type = "source")

# 4) CRAN-like check
rcmdcheck::rcmdcheck(
  tarball,
  args = c("--as-cran","--run-donttest","--timings"),
  error_on = "warning",
  check_dir = "."
)
'

rm -f gap/src/*.so gap/src/*.o || true
