#ifndef PGC_H
#define PGC_H

#define MAX_LOC 30

enum boolean{false,true};
#define _swap_(a,b) __swap__(a,b)
void __swap__(int *a,int *b) {int t;t=*a;*a=*b;*b=t;}
struct individual
{
  char id[20];
  int affection;
  int locus[MAX_LOC][2];
  int gtype[MAX_LOC];
} p_t;

static int n_obs=0,n_loci=0,handle_missing=0,with_id;
int getloci(char*),getdat(char*),getlocim(char*),getdatm(char*);
int a2g(int,int);
int g2a(int,int*,int*,int*);
int nloci,alleles[MAX_LOC],permute,npermute;
double nall[MAX_LOC],np[MAX_LOC],nnp[MAX_LOC],position(int,int*,int),positionm(int,int*,int);
int cc=false,sel[MAX_LOC],selp[MAX_LOC],isgenotype,iogenotype;
int selected,selectn,selectp,selidx[MAX_LOC],selndx[MAX_LOC],selpdx[MAX_LOC];
float freq,pen0,pen1,pen2;
int sample_size,cases;

typedef struct
{
  int l[MAX_LOC],u[MAX_LOC];
} phenotype;
phenotype *alist;

typedef struct node_type
{
  double genid;           /*unique genotype identifier*/
  int nca;                  /*number of cases*/
  int nco;                  /*number of controls*/
  int l[MAX_LOC],u[MAX_LOC];  /*actual phenotypes*/
  struct node_type *left, *right;
} node;

node *rt=0;
node *itree(node*,double);
node *stree(node*,double);
node *dtree(node*,double);
void inorder(node*),preorder(node*),postorder(node*);
void ctree(node*, double*, int*),rtree(node*),ptree(node*,int,FILE*);

#define maxloci MAX_LOC
#define maxalleles 50
#define M maxalleles*(maxalleles+1)/2

int digits=maxloci;
enum {CONTROL, CASE};

typedef struct newrec {
        int id,cc;
        int k[maxloci],locus[maxloci][2];
        struct newrec *next;
        } *list;
list Last, r;

list rsort(list, int), rsort1(list);
int getsize(FILE*), noid(char*,FILE*);

#endif
