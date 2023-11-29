#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void get_tdt_similarity(void *, void *, void *, void *);
extern void hap_read(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void hap_transmit(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void set_tdt_similarity(void *, void *, void *, void *);
extern void tdt_quad(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"get_tdt_similarity", (DL_FUNC) &get_tdt_similarity,  4},
    {"hap_read",           (DL_FUNC) &hap_read,           13},
    {"hap_transmit",       (DL_FUNC) &hap_transmit,       14},
    {"set_tdt_similarity", (DL_FUNC) &set_tdt_similarity,  4},
    {"tdt_quad",           (DL_FUNC) &tdt_quad,           10},
    {NULL, NULL, 0}
};

void R_init_tdthap(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
