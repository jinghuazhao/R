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
#define time timexxx
#include "S.h"
#undef time

/*
** Memory defined with S_alloc is removed automatically by S.
**  That with "CALLOC" I have to remove myself.  Use the
**  latter for objects that need to to persist between calls
*/
#if( defined(SPLUS_VERSION) && SPLUS_VERSION >= 5000)
#define ALLOC(a,b)  S_alloc(a,b,S_evaluator)
#define CALLOC(a,b) S_ok_calloc((size_t)(a), b, S_evaluator)
#else
#define ALLOC(a,b)  S_alloc(a,b)
#define CALLOC(a,b) S_ok_calloc((unsigned)(a), b)
#endif

/*
** This next is to make it easier to have common code with R, where
** "integers" are int.  In S they are long.
*/
#ifndef Sint
#define Sint int
#endif
