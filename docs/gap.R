#!/usr/bin/env bash
set -euo pipefail

cd ~/R
module load ceuadmin/R
export PATH="/home/jhz22/R-devel/bin:$PATH"
export MAKEFLAGS="-j$(nproc)"
export _R_CHECK_CRAN_INCOMING_=TRUE
export _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE
export _R_CHECK_FORCE_SUGGESTS_=FALSE
export _R_CHECK_LIMIT_CORES_=TRUE
export _R_CHECK_BUILD_VIGNETTES_=FALSE

Rscript -e '
cat(R.version.string, "\n")
pkg <- "gap"
devtools::document(pkg)
dir.create("gap/inst/doc", recursive=TRUE, showWarnings=FALSE)
rmarkdown::render("vignettes/gap.Rmd",
           output_file="gap_incl.html", output_dir="gap/inst/doc", quiet=TRUE)
tarball <- pkgbuild::build(pkg, clean=TRUE, args=c("--compact-vignettes=both"))
install.packages(tarball, repos=NULL, type="source")
rcmdcheck::rcmdcheck(
  tarball,
  args = c("--as-cran", "--run-donttest", "--timings"),
  error_on = "warning"
)
'

rm -f gap/src/*.{o,so} 2>/dev/null || true
