#ifndef BDSS_H
#define BDSS_H

#include <R.h>

/* Compatibility typedef retained for the original sources. */
#ifndef Sint
#define Sint int
#endif

/* Allocate memory that is automatically released when the
   .C/.Call invocation returns. */
#define ALLOC(a, b) R_alloc((a), sizeof(b))

#endif
