#ifndef NGHDS_H
#define NGHDS_H

#include "gdef.h"

typedef struct _Ind {
	struct _Ind *pa,*ma;
	int self;
	int index;
	int component;
	struct _Marlist *marriages;
	int ntrips;
	int *trips;
	int nbips;
	int *bips;
} Ind;

typedef struct _Component {
	int index;
	struct _Indlist *people;
} Component;

typedef struct _Mar {
	struct _Ind *spouse;
	struct _Indlist *kids;
} Mar;

typedef struct _Marlist {
	struct _Mar *current;
	struct _Marlist *next;
} Marlist;

typedef struct _Indlist {
	struct _Ind *current;
	struct _Indlist *next;
} Indlist;

#ifdef DEFINE_DATA
#define GLOBAL
#else
#define GLOBAL extern
#endif

#define INVALID_INDEX -1

#define UNKNOWN 0	/* Unknown pedigree names are 0 by convention*/

#define MALE   1	/* "Usual" coding for sex */
#define FEMALE 2

#define MAX_NAMES 100000	/* default for biggest name	*/
#define MAX_NCOMPS   100	/* initial value only for 	*/
				/* maximum no. of components	*/

GLOBAL int max_names;	/* was MAX_NAMES */
GLOBAL int **posit;	/* was int posit[MAX_NAMES][2] */
GLOBAL int ncol_data;	/* was 7 (or 8) */
GLOBAL int **data;	/* was data[MAX_SIZE][7] */
GLOBAL int pedsize;
GLOBAL int namestring;
GLOBAL int n_compnts;	/* # of components found	*/
GLOBAL int max_ncmps;	/* # of components allocated	*/
GLOBAL int max_cpsize;  /* maximum component size found */

GLOBAL Ind *nghd;	   
GLOBAL Component *cmplist;  /* vector of components */
GLOBAL Ind nullnode;

void alloc_ped(void);
void alloc_pos(void);
void free_ped(void);
void free_nghds(void);
void free_cmplists(void);
void position(int **, int **);
void nghds(int **);

int count_founders(void);
int input_ped(FILE *, int, int, int, int *, int *, int **, double **);
int intmax(int *aa, int n1, int n2);

double kinship(Ind *, Ind *), inbreeding(Ind *);

#endif
