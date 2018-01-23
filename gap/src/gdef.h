/*  contains general data and function definitions	*/

#ifndef GDEF_H
#define GDEF_H 

#include <ctype.h>
#include <stdarg.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

/* #include "dmalloc.h" */ 			/* for use with	dmalloc	*/ 
						/* debugger only	*/
#ifdef  DEFINE_DATA
#define GLOBAL 
#else
#define GLOBAL extern
#endif /* DEFINE_DATA */

#define EXIT_MALLOC 44			/*exit status code for malloc	*/
					/*failure			*/

GLOBAL time_t start_time;		/*start time for program executn*/

GLOBAL FILE * pfile;			/*pointer to current input file	*/

GLOBAL int gethostname (char *, int);
GLOBAL void give_up(char *, int, char *, ...);
GLOBAL void print_head (char *);
GLOBAL void print_tail ();
GLOBAL void print_val (double);
GLOBAL void print_val5 (double);
GLOBAL void alloc_chk (void *, char *);

#endif
