\name{jma.cojo}
\docType{data}
\alias{jma.cojo}
\title{A data containing independent GWAS hits as from GCTA}
\description{This data is used to illustrate cis/trans classification, containg the following columns:
\preformatted{
    prot Chr                SNP       bp refA       freq          b        se
1 4E.BP1  19 chr19:54327313_A_C 54327313    A 0.20550900  0.4510040 0.0243056
2 4E.BP1  19 chr19:54329063_G_T 54329063    T 0.10023500 -0.3233240 0.0333274
3    ADA  19 chr19:54327313_A_C 54327313    A 0.20550900  0.3542660 0.0246266
4    ADA  20 chr20:37456819_C_T 37456819    T 0.00388582 -0.2473080 0.1749800
5    ADA  20 chr20:38196991_G_T 38196991    G 0.00236927 -0.0171435 0.2238980
6    ADA  20 chr20:38603207_A_G 38603207    A 0.17074600 -0.0269075 0.0271976
            p       n  freq_geno         bJ     bJ_se           pJ        LD_r
1 2.48545e-74 6483.69 0.20079500   0.426476 0.0251676  2.07907e-64 -0.13397800
2 4.69307e-22 6480.60 0.08846920  -0.246444 0.0338712  3.44090e-13  0.00000000
3 5.47833e-46 6441.97 0.20079500   0.354266 0.0250171  1.59869e-45  0.00000000
4 1.57618e-01 5553.51 0.00497018  -5.873090 0.2241210 2.32892e-151 -0.00633091
5 9.38970e-01 5556.57 0.00198807 -13.473100 0.3790980 1.18609e-276  0.02467370
6 3.22550e-01 6285.16 0.15009900  -0.299797 0.0278787  5.69806e-27  0.11116200
  UniProtID
1    Q13541
2    Q13541
3    P00813
4    P00813
5    P00813
6    P00813
}
}
\usage{data(jma.cojo)}
\format{A data frame containing 445 rows and 16 columns}
\source{Undisclosed}
  
\keyword{datasets}
