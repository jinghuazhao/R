#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void g2a_(int *s, int *l, int *u, int *t);
extern void gcontrol_c(double *kdata, int *nkdata, double *zeta, double *kappa,
                       double *tau2, double *epsilon, int *ngib, int *burn,
                       int *idumR, double *deltot, double *x, double *A);
extern void gcx(int *verbose,int *Rhandlemissing,int *Rconvll,double *Reps,double *Rtol,
                int *Rmaxit,double *Rpl,double *precis,int *gid,int *Rnloci,int *Rloci,
                int *Robscom,int *Rhapall,int *genotype,int *count,int *Rxdata,int *sex,
                int *hapid,double *prob,double *Rh0,double *Rh1,double *lnl0,double *lnl1,
                int *npusr,int *npdat,int *iter,int *converge,char **assignment);
extern void gif_c(int *data, int *famsize, int *gifset, int *giflen, double *gifval);
extern void hap_c(int *nobs, char **idstr, char **gdata, int *nloci, int *loci,
                  char **names, double *mb, double *pr, double *po, double *to,
                  double *th, int *maxitt, int *n, int *sst, int *rst, int *rp,
                  int *ro, int *rv, double *sd, int *mm, int *mi, int *mc,
                  double *ds, double *de, int *q, double *l1, int *niter,
                  int *converged, char **hapfile, char **assignfile);
extern void hwe_hardy(int *a, int *alleles, int *seed, int *gss,
                      double *p, double *se, double *swp);
extern void kbyl(int *nalleles1, int *nalleles2, double *h, double *haplotypes,
                 double *Dp, double *VarDp, double *Dijtable, double *VarDijtable,
                 double *X2table,
                 double *Dmaxtable, double *Dijptable, double *VarDijptable,
                 double *x2, double *seX2, double *rho, double *seR, int *optrho,
                 double *klinfo, int *verbose);
extern void kbylem(double *obs, int *nalleles1, int *nalleles2,
                    double *Rh, double *l0, double *l1);
extern void kin_morgan(int *data, int *pedsize, int *pedindex, double *kin);
extern void makeped_c(char **pifile, char **pofile, int *autoselect,
                      int *withloop, char **loopfile,
                      int *autoproband, char **probandfile);
extern void mia_c(char **hapfile, char **assfile, char **miafile,
                  int *so, int *ns, int *mi, int *allsnps, int *sas);
extern void onelocus(float *y1, float *p1);
extern void pgc_c(int *gdata, int *handlemiss, int *nobs, int *nloci,
                  int *alleles, int *wt, int *gret, int *withid,
                  double *idsave, int *obscom);
extern void score_all(int *allele, int *nn, double *nplscore);
extern void score_pairs(int *allele, int *nn, double *nplscore);
extern void tbyt(double *h, double *haplotypes, double *D, double *VarD,
                 double *Dmax, double *VarDmax, double *Dprime, double *VarDprime,
                 double *x2, double *lor, double *vlor);
extern void twolocus(float *y12, float *p1, float *p2);
extern void x22k(int *a, int *tablen, double *x2a, double *x2b, int *col1, int *col2, double *p);

/* .Fortran calls */
extern void F77_NAME(family_)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(runifamily)(void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"g2a_",        (DL_FUNC) &g2a_,         4},
    {"gcontrol_c",  (DL_FUNC) &gcontrol_c,  12},
    {"gcx",         (DL_FUNC) &gcx,         28},
    {"gif_c",       (DL_FUNC) &gif_c,        5},
    {"hap_c",       (DL_FUNC) &hap_c,       30},
    {"hwe_hardy",   (DL_FUNC) &hwe_hardy,    7},
    {"kbyl",        (DL_FUNC) &kbyl,        19},
    {"kbylem",      (DL_FUNC) &kbylem,       6},
    {"kin_morgan",  (DL_FUNC) &kin_morgan,   4},
    {"makeped_c",   (DL_FUNC) &makeped_c,    7},
    {"mia_c",       (DL_FUNC) &mia_c,        8},
    {"onelocus",    (DL_FUNC) &onelocus,     2},
    {"pgc_c",       (DL_FUNC) &pgc_c,       10},
    {"score_all",   (DL_FUNC) &score_all,    3},
    {"score_pairs", (DL_FUNC) &score_pairs,  3},
    {"tbyt",        (DL_FUNC) &tbyt,        11},
    {"twolocus",    (DL_FUNC) &twolocus,     3},
    {"x22k",        (DL_FUNC) &x22k,         7},
    {NULL, NULL, 0}
};

static const R_FortranMethodDef FortranEntries[] = {
    {"family_",    (DL_FUNC) &F77_NAME(family_),    9},
    {"runifamily", (DL_FUNC) &F77_NAME(runifamily), 7},
    {NULL, NULL, 0}
};

void R_init_gap(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
