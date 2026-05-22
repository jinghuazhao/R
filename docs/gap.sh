#!/usr/bin/env bash
set -euo pipefail

cd ~/R
module load ceuadmin/R
export PATH="/home/jhz22/R-devel/bin:$PATH"
export _R_CHECK_CRAN_INCOMING_=TRUE
export _R_CHECK_CRAN_INCOMING_REMOTE_=TRUE
export _R_CHECK_FORCE_SUGGESTS_=FALSE
export _R_CHECK_LIMIT_CORES_=TRUE
export MAKEFLAGS="-j$(nproc)"

Rscript -e 'cat(R.version.string,"\n");devtools::document("gap")'
Rscript -e 'rmarkdown::render("vignettes/gap.Rmd",
            output_file="gap_incl.html", output_dir="gap/inst/doc", clean=TRUE, quiet=TRUE)'
R CMD build --compact-vignettes=both gap
TARBALL=$(ls -t gap_*.tar.gz | head -1)
R CMD INSTALL "$TARBALL"
R CMD check --as-cran --run-donttest "$TARBALL"

rm -f gap/src/*.{o,so} 2>/dev/null || true
