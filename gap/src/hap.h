#define VERSION "0.2.1"
#define MAX_LINE_LEN 1024
#define MAX_ID_LEN 64
#define MAX_FILENAME_LEN 128

typedef struct {
  char* id;
  double prior, posterior;
  short *loci;
} HAP;

typedef struct {
  int anum;
  char one[3], two[3];
} CODE;

/* Function declarations */

HAP* new_hap(char *id, double prior, double posterior, char *loci);
HAP* cpy_hap(HAP *old);
void kill_hap(HAP* old);
int cmp_hap(const void *one, const void *two);
int more_probable(const void *one, const void *two);
long hap_expand(long n_hap, long max_haps, HAP** so_list, int random);
void hap_prior(long n_hap, HAP** ho_list, double min);
void hap_prior_restart(long, HAP**);
void hap_prior_restore(long, HAP**, double *);
void sample_prior(long n_hap, HAP** ho_list, double prior_df);
long hap_posterior(long n_hap, HAP** so_list, double min, double *llh, int);
void hap_posterior_restart(long, HAP**);
void sample_posterior(long, HAP**);
void hap_list(FILE *, long, CODE *, HAP**);
long n_unique_haps(long n_hap, HAP **ho_list);
void unique_haps(long n_hap, HAP **ho_list, HAP **unique);
long check_hap(long n_hap, HAP **);
int encode(char *, CODE *);
int gt_read(int, char**, char**, int *, CODE *, HAP **, HAP **);
int allele_code(int, CODE);
void ranord(int, int*);
int talloc(long);
long memavail(int);
int hap_write(FILE *, int, char **, CODE *, int *, long, HAP **, int, double,
               int, int);
double rangamma(double);

/* Unix/PC random number functions, 08/05/2002 JH Zhao */

#ifndef SEED /*srand48*/
/*
  #define SEED srand
*/
#endif
#ifndef URAN /*drand48*/
/*
  #define URAN() rand()/(double)RAND_MAX
*/
  #define URAN() unif_rand()
#endif

/* JH Zhao 1/6/1999 IoP */

double **allocateU(int []);
void freeU(double**);
