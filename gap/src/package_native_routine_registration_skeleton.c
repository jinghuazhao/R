#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void g2a_(void *, void *, void *, void *);
extern void gcontrol_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gcx(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gif_c(void *, void *, void *, void *, void *);
extern void hap_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hwe_hardy(void *, void *, void *, void *, void *, void *, void *);
extern void kbyl(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void kbylem(void *, void *, void *, void *, void *, void *);
extern void kin_morgan(void *, void *, void *, void *);
extern void makeped_c(void *, void *, void *, void *, void *, void *, void *);
extern void mia_c(void *, void *, void *, void *, void *, void *, void *, void *);
extern void onelocus(void *, void *);
extern void pgc_c(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void score_all(void *, void *, void *);
extern void score_pairs(void *, void *, void *);
extern void tbyt(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void twolocus(void *, void *, void *);
extern void x22k(void *, void *, void *, void *, void *, void *, void *);

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

