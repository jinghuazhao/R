#!/usr/bin/env bash
set -euo pipefail

cd ~/R
module load ceuadmin/R

function cli()
# Reference workflow
{
version=$(awk '/^Version:/{print $2}' ~/R/gap/DESCRIPTION)

R CMD build --compact-vignettes=both --md5 --resave-data --log gap
R CMD INSTALL --compact-docs --data-compress=xz gap_${version}.tar.gz
R CMD check --as-cran --run-donttest gap_${version}.tar.gz

if command -v valgrind >/dev/null; then
  R CMD check --use-valgrind gap_${version}.tar.gz
fi
}

# CRAN env
export _R_CHECK_CRAN_INCOMING_=TRUE
export _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE
export _R_CHECK_FORCE_SUGGESTS_=FALSE
export _R_CHECK_LIMIT_CORES_=TRUE
export MAKEFLAGS="-j$(nproc)"

echo "Using R:"
Rscript -e 'cat(R.version.string,"\n")'

# Modern workflow
Rscript -e '
  devtools::document("gap")
  tarball <- pkgbuild::build("gap")
  install.packages(tarball, repos = NULL, type = "source")
  library(gap)
  rcmdcheck::rcmdcheck(
    tarball,
    args = c("--as-cran","--run-donttest"),
    error_on = "warning"
  )
'
rm -f gap/src/*.so gap/src/*.o || true
