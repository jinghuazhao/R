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
            output_file="gap.html", output_dir="vignettes", clean=TRUE, quiet=TRUE);
            knitr::purl("vignettes/gap.Rmd", output="vignettes/gap.R", documentation=0)'
R CMD build --compact-vignettes=both gap
find vignettes -maxdepth 1 -type f \( -name '*.png' -o -name '10081*' \) ! -name 'IL-12B_mhtplot.trunc.png' -delete
TARBALL=$(ls -t gap_*.tar.gz | head -1)
R CMD INSTALL "$TARBALL"
R CMD check --as-cran --run-donttest "$TARBALL"

rm -f gap/src/*.{o,so} 2>/dev/null || true
