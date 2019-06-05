# Read a GAS or LINKAGE format pedigree, return a digraph in the dot language 
# call dot to make pedigree drawing
# 
AWK=/bin/gawk
DOTEXE=/usr/local/bin/dot
# cygwin
# AWK=/bin/gawk
# DOTEXE=c:/local/graphviz/bin/dot

for fil in $*
do
  for ped in `$AWK '!/^[!#]/ {print $1}' $fil | sort -u`
  do
     echo "Pedigree $ped"
     $AWK -v ped=$ped '
     BEGIN { shape["m"]="box,regular=1"
             shape["1"]="box,regular=1"
             shape["f"]="circle"
             shape["2"]="circle"
             shade["y"]="blue"
             shade["2"]="blue"
             shade["n"]="grey"
             shade["1"]="grey"
             shade["x"]="green"
             shade["0"]="green"
     }
     !/^[!#]/ && $1==ped {
             sex[$2]=$5
             aff[$2]="x" ; if ($6 ~ /[012nyx]/) aff[$2]=$6
             if($3!="x" && $3!="0") {
               marriage[$3,$4]++
               child[$3,$4,marriage[$3,$4]]=$2
             }
     }
     END   { print "digraph Ped_" ped " {"
             print "# page =\"8.2677165,11.692913\" ;"
             print "ratio =\"auto\" ;"
             print "mincross = 2.0 ;"
             print "label=\"Pedigree " ped "\" ;"
             print "rotate=90 ;"
             for(s in sex) {
               print "\"" s "\" [shape=" shape[sex[s]] ","  \
                     " style=filled,color=" shade[aff[s]] "] ;"
             }
             for(m in marriage) {
               n=split(m,par,"\034")
               mating="\"" par[1] "x" par[2] "\""
               print mating "[shape=diamond,style=filled," \
                     "label=\"\",height=.1,width=.1] ;"
               print "\"" par[1] "\" -> " mating " [dir=none, weight=1] ;"
               print "\"" par[2] "\" -> " mating " [dir=none, weight=1] ;"
               for(k=1;k<=marriage[par[1],par[2]];k++) {
                 print  mating " -> \"" child[par[1],par[2],k] "\"" \
                        " [dir=none, weight=2] ;"
               }
             }
             print "}"
     }' $fil > $ped.dot 
     $DOTEXE -Tps $ped.dot -o $ped.ps
  done
done
$DOTEXE -Tps *.dot -o $*.ps
