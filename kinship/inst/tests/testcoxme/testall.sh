 cat setup.R > temp
 cat test0.R >> temp
 cat igchol.R stest0.R >> temp
 cat test1.R >> temp
 cat test2.R  >> temp
 cat test3.R  >> temp
 cat simple.R simple2.R >> temp
 cat frailty.kin.R >> temp
 cat ftest1.R >> temp
 cat ftest2.R >> temp
 cat ftest3.R >> temp
 cat bdstest.R >> temp
 cat readibd.R >> temp
 echo 'q()' >> temp

 R --save <temp >testall.out &

