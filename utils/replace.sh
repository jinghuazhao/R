#!/usr/bin/bash

# to replace heavily compressed files with default to enable FireFox

for f in jss.pdf gap.examples.pdf rnews.pdf pedtodot.pdf lmm-tr.pdf pan-tr.pdf tdthap-paper.pdf
do
  cp work/${f} vignettes
done

git add utils
git commit -m "utils"
git add vignettes
git commit -m "Package vignettes"
git push
