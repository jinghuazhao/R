#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(cebayes)(int *m, int *p, int *q, double *b, double *u, double *a, double *xtwxinv, double *sigma2, double *ztvinvx, int *g, double *wkgg, double *wkqp, double *wkpg, double *wkqg, double *wkp, double *varbeta, double *varb, double *covbbeta, double *wkpg2, double *wkqg2, double *xi, double *wkqq1, int *ntot, int *err);
extern void F77_NAME(ecmeml)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *iflag, int *err, int *msg, double *u, int *iter, int *sflag, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *cvgd, double *obeta, double *oxi, int *maxits, double *llvec, double *eps);
extern void F77_NAME(ecmerml)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *iflag, int *err, int *msg, double *u, int *iter, int *sflag, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *cvgd, double *obeta, double *oxi, int *maxits, double *llvec, double *eps, double *ztvinvx, double *a, double *wkqp);
extern void F77_NAME(fastmcmc)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *iflag, int *err, int *msg, double *u, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *maxits, double *abc, double *dinv, double *sqrtu, double *sigma2s, double *psis, int *g, double *wkgg, double *wkgg2, double *wkg, double *sig2hat, double *xihat, double *xigibbs, int *reject, double *ratios, double *df);
extern void F77_NAME(fastml)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *iflag, int *err, int *msg, double *u, int *iter, int *sflag, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *cvgd, double *oxi, int *maxits, double *llvec, double *eps, double *xiecme, int *g, int *reject, double *wkg, double *wkgg);
extern void F77_NAME(fastmode)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *isflags, int *err, int *msg, double *u, int *iter, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *cvgd, double *oxi, int *maxits, double *llvec, double *eps, double *xiecme, int *g, int *reject, double *ztvinvx, double *a, double *wkqp, double *wkg, double *wkgg, double *abc, double *dinv);
extern void F77_NAME(fastrml)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *iflag, int *err, int *msg, double *u, int *iter, int *sflag, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *cvgd, double *oxi, int *maxits, double *llvec, double *eps, double *xiecme, int *g, int *reject, double *ztvinvx, double *a, double *wkqp, double *wkg, double *wkgg);
extern void F77_NAME(mgibbs)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *w, double *vinv, int *pcol, double *pred, int *q, int *zcol, double *ztvinv, double *ztvinvz, int *isflag, int *err, int *msg, double *u, double *sigma2, int *p, int *xcol, double *beta, double *y, double *delta, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *wkqq1, double *wkqq2, double *xi, double *wkqnm, double *b, int *maxits, double *abc, double *dinv, double *sqrtu, double *sigma2s, double *psis);
extern float F77_NAME(rngs)(int *seed);

static const R_FortranMethodDef FortranEntries[] = {
    {"cebayes",  (DL_FUNC) &F77_NAME(cebayes),  24},
    {"ecmeml",   (DL_FUNC) &F77_NAME(ecmeml),   43},
    {"ecmerml",  (DL_FUNC) &F77_NAME(ecmerml),  46},
    {"fastmcmc", (DL_FUNC) &F77_NAME(fastmcmc), 50},
    {"fastml",   (DL_FUNC) &F77_NAME(fastml),   47},
    {"fastmode", (DL_FUNC) &F77_NAME(fastmode), 51},
    {"fastrml",  (DL_FUNC) &F77_NAME(fastrml),  50},
    {"mgibbs",   (DL_FUNC) &F77_NAME(mgibbs),   41},
    {"rngs",     (DL_FUNC) &F77_NAME(rngs),      1},
    {NULL, NULL, 0}
};

void R_init_lmm(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
