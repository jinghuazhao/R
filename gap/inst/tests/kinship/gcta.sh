#4-12-2013 MRC-Epid JHZ

gcta64 --grm 51 --pca 51 --out 51
gcta64 --reml --grm 51 --pheno 51.txt --out 51_qt
gcta64 --reml --grm 51 --pheno 51.txt --mpheno 2 --prevalence 0.25 --out 51_bt_25
gcta64 --reml --grm 51 --pheno 51.txt --mpheno 2 --out 51_bt


