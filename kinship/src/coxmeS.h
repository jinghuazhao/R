/* $Id: coxmeS.h,v 1.1 2002/12/30 21:53:43 therneau Exp $ */
/*
**   The S.h file defines a few things that I need, and hundreds that I don't.
** In particular, on some architectures, it defines a variable "time"
** which of course conflicts with lots of my C-code, 'time' being a natural
** variable name for survival models.
**   Thanks to Brian Ripley for suggesting a machine independent way of
** fixing this.
**
** The S_alloc function changed it's argument list from version 4 to 5, and
**   the ALLOC macro allows me to have common C code for the two versions,
**   with only this file "survS.h" changed.
*/
#ifndef COXMES_H
#define COXMES_H

#define time timexxx
#include <R.h>
#include <R_ext/Memory.h>
#undef time

/*
** Memory allocated with R_alloc() is automatically released when the
** .C/.Call invocation returns. Memory allocated with Calloc() must be
** released explicitly with Free().
*/
#define ALLOC(a, b) ((b *) R_alloc((R_xlen_t)(a), sizeof(b)))

#ifndef Calloc
#define Calloc(n, type) ((type *) R_Calloc((R_SIZE_T)(n), type))
#endif

#ifndef Free
#define Free(p) R_Free(p)
#endif

/*
** Compatibility typedef retained for the original sources.
*/
#ifndef Sint
#define Sint int
#endif

#endif
