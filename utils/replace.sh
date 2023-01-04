#!/usr/bin/bash

# to replace heavily compressed files with default to enable FireFox

for f in jss.pdf # gap.examples.pdf rnews.pdf
do
  cp work/${f} vignettes
done

git add vignettes
git commit -m "Package vignettes"
git push
