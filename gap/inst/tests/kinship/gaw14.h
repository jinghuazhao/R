/*16-05-2004 extracted from gaw14.c to avoid code duplication*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define MAX_FAM 150
#define MAX_OBS 1700
#define MAX_LOC 50
#define MAX_ALL 40
#define LINELEN 3500
#define DEBUG
  
typedef struct {int fid,id,allele[MAX_LOC][2];} LOC;
typedef struct {char name[10];int n,size[MAX_ALL];float freq[MAX_ALL];} mark;
typedef struct {int pid,id,fid,mid; char sex[2]; int aff;} ind;/**/
LOC *genotype;
mark locus[MAX_LOC];
ind *person;
int fsize[MAX_FAM],sample_size,families,nloci;
int clips[27]={125,127,129,140,147,149,151,
               434,439,440,441,442,
               718,719,725,726,727,728,729,730,731,732,733,734,735,743,744};
float map[MAX_LOC];
FILE *par,*ped,*gh,*makeped,*lop;

int getphen(),getfrq(char *),getmap(char *),getdat(char *),a2a();
#define haldane(x) (0.5*(1-exp(-2*x)))
