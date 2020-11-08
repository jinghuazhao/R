#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(bdiag)(int *q, int *m, double *sig);
extern void F77_NAME(ecme3)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vi, double *vh, int *pcol, double *pred, int *q, int *zcol, double *ztv, double *sig0, int *iflag, double *sig, double *psi, double *sigma2, int *p, int *xcol, double *beta, double *wkq1, double *wkq2, double *wkq3, double *y, double *delta, double *b, double *wk, double *w, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *llk, double *vmax, int *sflag, double *eps, double *obeta, double *opsi, int *maxits, int *iter, int *cvgd);
extern void F77_NAME(mgibbs)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *pcol, double *pred, int *q, int *zcol, double *ztz, int *patt, int *nstar, int *r, double *y, int *p, int *xcol, int *npatt, int *rmat, int *sflag, double *beta, double *sigma, double *psi, double *b, double *xtxinv, double *wkpp, double *wkpr, double *eps, double *wkrqrq1, double *wkrqrq2, double *sig, double *wkrr1, double *wkrr2, int *iter, double *wkqr1, double *wkqr2, double *wkqrv, int *nhyp, double *hyp, double *delta, int *iposn, int *pstfin, double *betas, double *sigmas, double *psis);
extern void F77_NAME(mgibbsbd)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *pcol, double *pred, int *q, int *zcol, double *ztz, int *patt, int *nstar, int *r, double *y, int *p, int *xcol, int *npatt, int *rmat, int *sflag, double *beta, double *sigma, double *psi, double *b, double *xtxinv, double *wkpp, double *wkpr, double *eps, double *wkrqrq1, double *wkrqrq2, double *sig, double *wkrr1, double *wkrr2, int *iter, double *wkqr1, double *wkqr2, double *wkqrv, int *nhyp, double *hyp, double *delta, int *iposn, int *pstfin, double *betas, double *sigmas, double *psis, double *wkqq1, double *wkqq2);
extern void F77_NAME(nopsi)(int *ntot, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vi, double *vh, int *pcol, double *pred, int *q, int *zcol, double *ztv, double *sig0, int *iflag, double *sig, double *psi, double *sigma2, int *p, int *xcol, double *beta, double *wkq1, double *wkq2, double *wkq3, double *y, double *delta, double *b, double *wk, double *w, double *xtw, double *xtwx, double *xtwy, double *xtwxinv, double *ll);
extern void F77_NAME(prelim)(int *ntot, int *subj, int *m, int *ist, int *ifin, int *occ, int *nmax, double *vmax, double *vh, double *vi, int *pcol, double *pred, int *q, int *zcol, double *ztv, double *sig0, int *iflag);
extern float F77_NAME(rngs)(int *seed);

static const R_FortranMethodDef FortranEntries[] = {
    {"bdiag",    (DL_FUNC) &F77_NAME(bdiag),     3},
    {"ecme3",    (DL_FUNC) &F77_NAME(ecme3),    43},
    {"mgibbs",   (DL_FUNC) &F77_NAME(mgibbs),   44},
    {"mgibbsbd", (DL_FUNC) &F77_NAME(mgibbsbd), 46},
    {"nopsi",    (DL_FUNC) &F77_NAME(nopsi),    34},
    {"prelim",   (DL_FUNC) &F77_NAME(prelim),   17},
    {"rngs",     (DL_FUNC) &F77_NAME(rngs),      1},
    {NULL, NULL, 0}
};

void R_init_pan(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
