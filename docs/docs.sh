# 14-3-2024 JHZ

cp -p ~/R/gap/ChangeLog ~/R/vignettes/ChangeLog.txt
for f in gap.html jss.pdf shinygap.html; do cp -p ~/hpc-work/R/gap/doc/${f} ~/R/vignettes; done
cp -p ~/hpc-work/R/pQTLdata/doc/pQTLdata.html ~/R/vignettes/
cp -p ~/R/gap.Rcheck/gap-manual.pdf ~/R/vignettes

git add .github
git commit -m ".github"
git add .gitignore
git commit -m "These are ignored"
git add CGR CGR_1.0-5.tar.gz
git commit -m "Classic Genetics in R"
git add kinship_1.1.4.tar.gz kinship_1.1.4.zip
git commit -m "kinship source (.tar.gz) and Windows (.zip) packages"
git add gap gap.enl
git commit -m "genetic analysis package"
git add gap.datasets
git commit -m "Datasets for 'gap'"
git add gap.examples
git commit -m "Examples for 'gap'"
git add lmm
git commit -m "linear mixed model"
git add pan
git commit -m "Multiple imputation for multivariate panel or clustered dat"
git add kinship
git commit -m "mixed-effects Cox models, sparse matrices, and modeling data from large pedigrees"
git add tdthap
git commit -m "Transmission/disequilibrium tests for extended marker haplotypes"
git add README.md
git commit -m "README"
git add docs
git commit -m "docs"
git add vignettes
git commit -m "vignettes"
git push
