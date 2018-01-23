#include <R.h>
#include "nghds.h"
#include <R.h>

#define max_size 1001
#define STRICT_CHECK

static void nullify(Ind *);
double kinship(Ind *, Ind *),inbreeding(Ind *);

Ind nullnode;

void kin_morgan(int *data, int *pedsize, int *pedindex, double *kin)
{
Ind *ind,*t1,*t2;
int i0,i1,i2,i,j,k,id,pa,ma;

nullify(&nullnode);
ind=(Ind *)malloc((max_size)*sizeof(Ind));
if (!ind)
{
  Rprintf("\nError to allocate memory for pedigree\n");
  return;
}
for (i=0;i<max_size;i++) nullify(&ind[i]);

for (j=1;j<=*pedsize;j++)
{
    i=j-1;
    id=data[i*3];
    pa=data[i*3+1];
    ma=data[i*3+2];
    i0=pedindex[i*3];
    i1=pedindex[i*3+1];
    i2=pedindex[i*3+2];
    if (i0) {ind[i0].self=id;ind[i0].index=i0;}
    if (i1) {ind[i1].self=pa;ind[i1].index=i1;}
    if (i2) {ind[i2].self=ma;ind[i2].index=i2;}
    ind[j].self=id;
    ind[j].index=i0;
}
for (j=1;j<=*pedsize;j++)
{
    i=j-1;
    pa=data[i*3+1];
    ma=data[i*3+2];
    i1=pedindex[i*3+1];
    i2=pedindex[i*3+2];
    t1=t2=&nullnode;
    if (pa) t1=&ind[i1];
    if (ma) t2=&ind[i2];
    ind[j].pa=t1;
    ind[j].ma=t2;
#ifdef STRICT_CHECK
    t1=&ind[i1];
    t2=&ind[i2];
    if ((pa && t1->self==UNKNOWN)||(ma && t2->self==UNKNOWN))
    {
       Rprintf("\nParents not in datafile, quit\n");
       Rprintf("pa=%5d ma=%5d t1->self=%5d t2->self=%5d\n",pa,ma,t1->self,t2->self);
       return;
    }
#endif
}
for (i=1;i<=*pedsize;i++)
{
    t2=&ind[i];
/*    Rprintf("%5d ",t2->self); */
}
/* Rprintf("\n"); */
k=0;
for (i=1;i<=*pedsize;i++)
{
    t1=&ind[i];
/*    Rprintf("%5d ",t1->self); */
    for (j=1;j<=i;j++)
    {
        kin[k]=kinship(&(ind[i]),&(ind[j]));
/*        Rprintf(" %f",kin[k]); */
        k++;
    }
/*    Rprintf("\n"); */
}

for (i=0;i<=*pedsize;i++) nullify(&ind[i]);
free(ind);

return;

}

static void nullify(Ind *nul)
/*
 * Function to make the nullnode for the neighbourhood routine.
 * Pointers to pa, ma and marriages are set to NULL. It is assumed that a
 * structure  of type Ind has already been set up as nullnode and defined
 * to be a global variable.
 */
{
        nul->self = UNKNOWN;
        nul->index = INVALID_INDEX;
        nul->pa = NULL;
        nul->ma = NULL;
        nul->marriages = NULL;
}

double kinship(Ind * a, Ind * b)
/*
 * Recursive program which returns the kinship coefficient between
 * individuals a and b.
 */
{
    if (a == &nullnode || b == &nullnode)
         return 0.0;

    if (a == b)
        return 0.5 + 0.5 * inbreeding(a);
    else
    if (a->pa->self == 0)
    {
        if (b->index < a->index)
            return 0.0;
        else
        if (b->pa->self == 0)
            return 0.0;
        else
            return (kinship(a, b->pa) + kinship(a, b->ma)) * 0.5;
    } else
    if (b->pa->self == 0)
    {
        if (a->index < b->index)
            return 0.0;
        else
            return (kinship(b, a->pa) + kinship(b, a->ma)) * 0.5;
    } else
    if (a->index < b->index)
        return (kinship(a, b->pa) + kinship(a, b->ma)) * 0.5;
    else
        return (kinship(b, a->pa) + kinship(b, a->ma)) * 0.5;
}

double inbreeding(Ind * a)
/*
 * A recursive program that returns the inbreeding coefficient
 * for individual a.
 */
{
    if (a == &nullnode)
         return 0.0;
    else
        return kinship(a->pa, a->ma);
}
