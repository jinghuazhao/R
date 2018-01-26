#ifndef bdsS_H
#define bdsS_H

#include "S.h"

/* The first line below is the correct line for Splus 
** the second the correct one for R */
/* #define Sint long*/
#define Sint int 

/* how to allocate memory */
#if( defined(SPLUS_VERSION) && SPLUS_VERSION >= 5000)
#define ALLOC(a,b) S_alloc(a,b, S_evaluator)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#endif

#endif
