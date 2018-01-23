/* 15-12-2015 MRC-Epid JHZ */

args s
use join.dta
set seed `s'
sample 70
sort id1
outsheet id1 id using ddjack.id, noname noquote replace
