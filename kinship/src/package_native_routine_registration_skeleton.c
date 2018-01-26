#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/*
  The following symbols/expressions for .NAME have been omitted

    routines[1]
    routines[2]
    ofile
    routines[3]

  Most likely possible values need to be added below.
*/

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void bdsmatrix_index1(void *, void *, void *, void *, void *, void *, void *, void *);
extern void bdsmatrix_index2(void *, void *, void *, void *);
extern void bdsmatrix_index3(void *, void *, void *);
extern void bdsmatrix_prod(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void bdsmatrix_prod3(void *, void *, void *, void *, void *, void *, void *);
extern void gchol(void *, void *, void *);
extern void gchol_bds(void *, void *, void *, void *, void *, void *);
extern void gchol_bdsinv(void *, void *, void *, void *, void *, void *, void *, void *);
extern void gchol_bdssolve(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void gchol_inv(void *, void *, void *);
extern void gchol_solve(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"bdsmatrix_index1", (DL_FUNC) &bdsmatrix_index1, 8},
    {"bdsmatrix_index2", (DL_FUNC) &bdsmatrix_index2, 4},
    {"bdsmatrix_index3", (DL_FUNC) &bdsmatrix_index3, 3},
    {"bdsmatrix_prod",   (DL_FUNC) &bdsmatrix_prod,   9},
    {"bdsmatrix_prod3",  (DL_FUNC) &bdsmatrix_prod3,  7},
    {"gchol",            (DL_FUNC) &gchol,            3},
    {"gchol_bds",        (DL_FUNC) &gchol_bds,        6},
    {"gchol_bdsinv",     (DL_FUNC) &gchol_bdsinv,     8},
    {"gchol_bdssolve",   (DL_FUNC) &gchol_bdssolve,   9},
    {"gchol_inv",        (DL_FUNC) &gchol_inv,        3},
    {"gchol_solve",      (DL_FUNC) &gchol_solve,      4},
    {NULL, NULL, 0}
};

void R_init_kinship(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
