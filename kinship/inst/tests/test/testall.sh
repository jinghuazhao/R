#
# Do the full set of tests
#
cat setup.R > temp
cat kin1.R >> temp
cat makefam.R >> temp
cat bdstest.R >> temp
cat chtest.R >> temp
cat chtest2.R >> temp
cat gtest.R   >> temp
cat matrix.R >>  temp
cat tinv.R   >>  temp
echo 'q()'   >> temp

R --save <temp >testall.out
