/*------------------------------------------------*/
/* calculation of kinship/inbreeding coefficients */
/* clone of Morgan v2.0                           */
/* JH Zhao 13-16/07/1999 IoP                      */
/* - 13/07 it works                               */
/* - 14/07 add nullify                            */
/* - 15/07 separate kin routine from main()       */
/* - 16/07 add temporary array to keep parents    */
/* - 20/07 add check to precedence of parents     */
/*------------------------------------------------*/

#include "nghds.h"

/*
 * maximum family size
 */
#define max_size 200
#define STRICT_CHECK

int kin(FILE*,FILE*);
static void nullify(Ind *);
double kinship(Ind *, Ind *),inbreeding(Ind *);

Ind nullnode;

int main(int argc, char *argv[])
{
FILE *fp,*fo;
int rt;

fp=stdin;
if (argc>1) fp=fopen(argv[1],"r");
if (!fp) fprintf(stderr,"\nNo such file called %s\n",argv[1]);

fo=stdout;
if (argc>2) fo=fopen(argv[2],"w");
if (!fo) fprintf(stderr,"\nCannot open output file %s\n",argv[2]);

rt=kin(fp,fo);

fclose(fp);
fclose(fo);

if (rt) return 1;
return 0;

}

int kin(FILE *fp, FILE *fo)
/* 
 * Routine to perform the calculation based on fp and output fo
 */
{
Ind *ped,*t1,*t2;
int pedsize,i,j,pid,id,pa,ma,items,**data;
char buffer[255];

data=(int **)malloc((max_size+1)*sizeof(int*));
if (!data)
{
  fprintf(stderr,"\nError to allocate memory for temporary data\n");
  return 1;
}
for (i=0;i<=max_size;i++)
{
    data[i]=(int *)malloc(5*sizeof(int));
    if (!data[i])
    {
       fprintf(stderr,"\nError to allocate memory for temporary arrays\n");
       return 1;
    }
}
nullify(&nullnode);
ped=(Ind *)malloc((max_size+1)*sizeof(Ind));
if (!ped)
{
  fprintf(stderr,"\nError to allocate memory for pedigree\n");
  return 1;
}
for (i=0;i<=max_size;i++)
{
     nullify(&ped[i]);
     for (j=0;j<5;j++)
         data[i][j]=0;
}

fprintf(fo,"\nThe original family (PID ID PA MA): \n\n");

i=0;
while (fgets(buffer,254,fp))
{
      items=sscanf(buffer,"%d %d %d %d",&pid,&id,&pa,&ma);
      if (items==3)
      {
         fprintf(stderr,"\nProbably you missed pedigree id\n");
         return 1;
      }
      fprintf(fo,"%4d %5d %5d %5d\n",pid,id,pa,ma);
      i++;
      data[i][1]=pid;
      data[i][2]=id;
      data[i][3]=pa;
      data[i][4]=ma;
      ped[i].self=id;
      ped[i].index=i;
#ifdef STRICT_CHECK
      t1=&ped[pa];
      t2=&ped[ma];
      if ((pa && t1->self==UNKNOWN)||(ma && t2->self==UNKNOWN))
      {
         fprintf(stderr,"\nParents not in datafile, quit\n");
         return 1;
      }
#endif
}
pedsize=i;
for (i=1;i<=pedsize;i++)
{
    t1=t2=&nullnode;
    pa=data[i][3];
    ma=data[i][4];
    if (pa) t1=&ped[pa];
    if (ma) t2=&ped[ma];
    ped[i].pa=t1;
    ped[i].ma=t2;
}
fprintf(fo,"\nThe kinship coefficient matrix:\n\n");
for (i=1;i<=pedsize;i++)
{
    fprintf(fo,"%5d ",i);
    for (j=1;j<=i;j++)
        fprintf(fo," %lf",kinship(&(ped[i]),&(ped[j])));
    fprintf(fo,"\n");
}

for (i=0;i<=max_size;i++) free(data[i]);
free(data);
for (i=0;i<=pedsize;i++) nullify(&ped[i]);
free(ped);

return 0;

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
