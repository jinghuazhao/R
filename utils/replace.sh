#!/usr/bin/bash

# to replace heavily compressed files with default to enable FireFox

for f in jss.pdf gap.examples.pdf rnews.pdf
do
  cp work/${f} vignettes
done

git add README.md
git commit -m "README"
git add utils
git commit -m "utils"
git add vignettes
git commit -m "Package vignettes"
git push
