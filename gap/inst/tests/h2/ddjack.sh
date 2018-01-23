#! 15-12-2015 MRC-Epid JHZ

for i in `seq 1 10000`
do
   stata -b do ddjack.do ${i}
   gcta64 --reml --reml-est-fix --grm chips/join --thread-num 20 --out ddjack/${i} --prevalence 0.05 \
          --pheno join.txt --gxe gxe.txt --keep ddjack.id  --qcovar qcovar.txt
done

