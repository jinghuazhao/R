#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* .C calls */
extern void get_tdt_similarity(int *, double *, int *, double *);
extern void hap_read(int *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, int *, char **);
extern void hap_transmit(int *, int *, int *, int *, int *, int *, int *, int *, double *, int *, int *, int *, int *, char **);
extern void set_tdt_similarity(int *, double *, int *, double *);
extern int tdt_quad(int *, int *, char **, int *, int *, int *, int *, int *, int *, double *);

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
