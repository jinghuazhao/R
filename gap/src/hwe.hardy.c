#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <Rmath.h>

#define  MAX_ALLELE    100
#define  LENGTH        MAX_ALLELE * ( MAX_ALLELE + 1 ) / 2
#define  STR_END       '\0'
#define  MIN(x, y)     ((x) < (y)) ? (x) : (y)
#define  RATIO(u, v)   ((double) (u) ) / ( 1.0 + (double) (v) )
#define  TRANS(x)      (MIN(1.0, x))/2.0
#define  LL(a, b)      a * ( a + 1 ) / 2  + b
#define  L(a, b)       ( a < b ) ? b*(b + 1)/2 + a : a*(a+1)/2 + b
/*
#define  drand48()     rand()/(double)RAND_MAX
*/
typedef struct {int i1, i2, j1, j2, type; double cst;} Index;
typedef struct {double p_value, se; int swch_count[3];} outcome;
typedef struct {int group, size, step;} randomization;

double cal_probn (int *, Index, double, int *);
double ln_p_value (int *, double);
double log_factorial (int);
double cal_const (int *n, int);
void cal_n (int *, int *);
void select_index (Index *, int);
void random_choose (int *, int *, int);
void test_switch (int *, Index, int *, int *, double *, double *);
void ndo_switch (int *, Index, int);
void stamp_time (time_t, FILE *);

int *work, no_allele;

void hwe_hardy(int *a, int *alleles, int *seed, int *gss,
         double *p, double *se, double *swp)
{
  int b[LENGTH], n[MAX_ALLELE];
  double ln_p_observed, ln_p_simulated, p_mean, p_square;
  double constant, p_simulated, total_step;
  int total, counter, actual_switch;
  Index index;
  randomization sample;
  outcome result={0,0,{0,0,0}};
  int i, j, k;

  GetRNGstate();
  sample.group = gss[0];
  sample.size = gss[1];
  sample.step = gss[2];
/*
  srand(*seed);
*/
  no_allele=*alleles;
  total = 0;
  for ( i = 0; i < no_allele; ++i ) 
    for ( j = 0; j <= i; ++j ) {
         k = LL(i, j);
         b[k] = a[k];
         total += b[k];
    }
  cal_n (b, n );
  constant = cal_const (n, total );
  ln_p_observed = ln_p_value ( b, constant );
  ln_p_simulated = ln_p_observed = 0.0;
  work=(int*)malloc(MAX_ALLELE*sizeof(int));
  if(!work) {
    REprintf("not enough memory\n");
    return;
  }
  for ( i = 0; i < sample.step; ++i ) {
    select_index ( &index, no_allele );
    ln_p_simulated = cal_probn(b, index, ln_p_simulated, &actual_switch);
    ++result.swch_count[actual_switch];
  }
  p_mean = p_square = 0.0;
  for ( i = 0; i < sample.group; ++i ) {
    counter = 0;
    for ( j = 0; j < sample.size; ++j ) {
      select_index ( &index, no_allele );
      ln_p_simulated = cal_probn(b, index, ln_p_simulated, &actual_switch);
      if ( ln_p_simulated <= ln_p_observed ) ++counter;
      ++result.swch_count[actual_switch];
    }
    p_simulated = (double) counter  / sample.size;
    p_mean += p_simulated;
    p_square += p_simulated * p_simulated;
  }
  free(work);
  p_mean /= sample.group;
  result.p_value = p_mean;
  result.se = p_square / (sample.group)/(sample.group - 1.0)
            - p_mean / (sample.group - 1.0) * p_mean;
  result.se = sqrt ( result.se );
  total_step = sample.step + sample.group * sample.size;
  *p=result.p_value;
  *se=result.se;
  swp[0]=result.swch_count[1]/total_step * 100;
  swp[1]=result.swch_count[2]/total_step * 100;
  swp[2]=(result.swch_count[1]+result.swch_count[2])/total_step * 100;
  PutRNGstate();
}

double cal_probn (int *a, Index index, double ln_p_old, int *actual_switch)
{
  double rand_num, p1_ratio=0, p2_ratio=0;
  double ln_p_new;
  int switch_ind, type=0;
  int k11, k22, k12, k21;

  k11 = L(index.i1, index.j1);
  k12 = L(index.i1, index.j2);
  k21 = L(index.i2, index.j1);
  k22 = L(index.i2, index.j2);

  *actual_switch = 0;
  switch_ind = 0;

  if ( index.type <= 1 ) {
    if ( a[k11] > 0 && a[k22] > 0 ) {
         switch_ind = 1;
         type = 0;
         p1_ratio = RATIO(a[k11], a[k12]) *  RATIO(a[k22], a[k21]) * index.cst;
    }
    if ( a[k12] > 0 && a[k21] > 0 ) {
         switch_ind += 1;
         type = 1;
         p2_ratio = RATIO(a[k12], a[k11]) *  RATIO(a[k21], a[k22]) / index.cst;
    }
  } else {
    if ( a[k11] > 0 && a[k22] > 0 ) {
         switch_ind = 1;
         type = 0;
         p1_ratio = RATIO(a[k11],a[k12] + 1.0)*RATIO(a[k22],a[k12]) * index.cst;
    }
    if ( a[k12] > 1 ) {
         switch_ind += 1;
         type = 1;
         p2_ratio = RATIO(a[k12],a[k11]) * RATIO(a[k12] - 1,a[k22]) / index.cst;
    }
  }
  switch (switch_ind)
    {
    case 0:
      ln_p_new = ln_p_old;
      break;
    case 1:
      if ( type == 1 ) p1_ratio = p2_ratio;
      rand_num = unif_rand();
      if ( rand_num < TRANS( p1_ratio ) ) {
           if ( type == 0 ) {
             --a[k11];
             --a[k22];
             ++a[k12];
             ++a[k21];
           } else {
             ++a[k11];
             ++a[k22];
             --a[k12];
             --a[k21];
           }
           ln_p_new = ln_p_old + log (p1_ratio);
           *actual_switch = 1;
         } else ln_p_new = ln_p_old;
       break;
    default:
      rand_num = unif_rand();
      if ( rand_num <= TRANS(p1_ratio)) {
        --a[k11];
        --a[k22];
        ++a[k12];
        ++a[k21];
        ln_p_new = ln_p_old + log (p1_ratio);
        *actual_switch = 2;
      } else if ( rand_num <= TRANS(p1_ratio) + TRANS(p2_ratio) ) {
        ++a[k11];
        ++a[k22];
        --a[k12];
        --a[k21];
        ln_p_new = ln_p_old + log (p2_ratio);
        *actual_switch = 2;
      } else ln_p_new = ln_p_old;
      break;
    }
  return (ln_p_new);
}

void select_index (Index *index, int no_allele)
{
  int i1, i2, j1, j2, k = 0, l = 0;

  random_choose ( &i1, &i2, no_allele );
  index->i1 = i1;
  index->i2 = i2;
  random_choose ( &j1, &j2, no_allele );
  index->j1 = j1;
  index->j2 = j2;
  if ( i1 == j1 ) ++k;
  if ( i1 == j2 ) ++k;
  if ( i2 == j1 ) ++k;
  if ( i2 == j2 ) ++k;
  index->type = k;
  if ( ( i1 == j1 ) || ( i2 == j2 ) ) ++l;
  index->cst = ( l == 1 ) ? pow(2.0, (double) k) : pow(2.0, - (double) k);
}

void random_choose (int *k1, int *k2, int k)
{
  int temp, i, not_find;

  for ( i = 0; i < k; ++i ) work[i] = i;
   *k1 = unif_rand() * k;
   --k;
  for ( i = *k1; i < k; ++i ) work[i] = i + 1;
  not_find = 1;
  while ( not_find ) {
    i = unif_rand() * k;
    *k2 = work[i];
    not_find = 0;
  }
  if ( *k1 > *k2 ) {
    temp = *k1;
    *k1 = *k2;
    *k2 = temp;
  }
}

void ndo_switch (int *a, Index index, int type)
{
  int k11, k22, k12, k21;

  k11 = L(index.i1, index.j1);
  k12 = L(index.i1, index.j2);
  k21 = L(index.i2, index.j1);
  k22 = L(index.i2, index.j2);

  if ( type == 0 ) {
    --a[k11];
    --a[k22];
    ++a[k12];
    ++a[k21];
  } else {
    ++a[k11];
    ++a[k22];
    --a[k12];
    --a[k21];
  }
}

double ln_p_value (int *a, double constant)
{
  register int i, j, l, temp;
  register double ln_prob;

  ln_prob = constant;
  temp = 0;
  for ( i = 0; i < no_allele; ++i ) {
    for ( j = 0; j < i; ++j ) {
         l = LL(i, j);
         temp += a[l];
      ln_prob -= log_factorial ( a[l] );
    }
         l = LL(i, i);
         ln_prob -= log_factorial ( a[l] );
  }
  ln_prob += temp * log ( 2.0 );
  return ( ln_prob );
}

void cal_n (int *a, int *n)
{
  register int i, j, l;

  for ( i = 0; i < no_allele; ++i ) {
    l = LL(i, i);
    n[i] = a[l];
    for ( j = 0; j < no_allele; ++j ) {
         l = L(i, j);
         n[i] += a[l];
    }
  }
}

double cal_const (int *n, int total)
{
  double constant;
  register int i;

  constant = log_factorial ( total ) - log_factorial ( 2*total );
  for ( i = 0; i < no_allele; ++i )
    constant += log_factorial ( n[i] );
  return ( constant );
}

double log_factorial (int k)
{
  register double result;

  if ( k == 0 ) result = 0.0;
  else result = log( (double)k ) + log_factorial ( k - 1 );
  return (result);
}

void stamp_time (time_t t1, FILE *outfile)
{
  time_t t2, now;

  time(&t2);
  fprintf (outfile, "\nTotal elapsed time: %.0f''\n", difftime(t2,t1));
  time(&now);
  fprintf (outfile, "Date and time: %s\n", ctime(&now));
}

#ifdef executable

int check_file (int, char *[], FILE **, FILE **);
int read_data (int *, int *, int *, randomization *, FILE **);
void print_data (int *, int, randomization, FILE **);
double cal_prob (int *, Index, double, int *);

int main(int argc, char *argv[])
{
  int a[LENGTH], n[MAX_ALLELE];
  double ln_p_observed, ln_p_simulated, p_mean, p_square;
  double constant, p_simulated, total_step;
  int no_allele, total, counter, actual_switch;
  Index index;
  randomization sample;
  outcome result;
  FILE *infile, *outfile;
  time_t t1;
  register int i, j;

  if ( check_file ( argc, argv, &infile, &outfile ) ) error("%d",1);
  time(&t1);
/*
  srand(3000);
*/
  if ( read_data ( a, &no_allele, &total, &sample, &infile ) ) error("%d",2);
  print_data ( a, no_allele, sample, &outfile );
  ln_p_observed = 0.0;
  ln_p_simulated = ln_p_observed;
  p_mean = p_square = (double) 0.0;
  result.p_value = result.se = (double) 0.0;
  result.swch_count[0] = result.swch_count[1] = result.swch_count[2] = 0;
  work=(int*)malloc(MAX_ALLELE*sizeof(int));
  if(!work) {
    REprintf("not enough memory\n");
    return;
  }
  for ( i = 0; i < sample.step; ++i ) {
    select_index ( &index, no_allele );
    ln_p_simulated = cal_prob(a, index, ln_p_simulated, &actual_switch);
    ++result.swch_count[actual_switch];
  }
  for ( i = 0; i < sample.group; ++i ) {
    counter = 0;
    for ( j = 0; j < sample.size; ++j ) {
    select_index ( &index, no_allele );
    ln_p_simulated = cal_prob(a, index, ln_p_simulated, &actual_switch);
    if ( ln_p_simulated <= ln_p_observed ) ++counter;
    ++result.swch_count[actual_switch];
  }
    p_simulated = (double) counter  / sample.size;
    p_mean += p_simulated;
    p_square += p_simulated * p_simulated;
  }
  p_mean /= sample.group;
  result.p_value = p_mean;
  result.se = p_square / ((double) sample.group)/(sample.group - 1.0)
    - p_mean / ( sample.group - 1.0 ) * p_mean;
  result.se = sqrt ( result.se );
  total_step = sample.step + sample.group * sample.size;
  fprintf(outfile, "Randomization test P-value: %7.4g  (%7.4g) \n",
                result.p_value, result.se);
  fprintf(outfile, "Percentage of partial switches: %6.2f \n",
                result.swch_count[1] / total_step * 100);
  fprintf(outfile, "Percentage of full switches: %6.2f \n",
                result.swch_count[2] / total_step * 100);
  fprintf(outfile, "Percentage of all switches: %6.2f \n",
                (result.swch_count[1] + result.swch_count[2]) / total_step * 100 );
  stamp_time ( t1, outfile );
  free(work);
}

int check_file (int argc, char *argv[], FILE **infile, FILE **outfile)
{

  int exit_value = 0;

  if ( argc != 3 ) {
    REprintf ("Bad commond.\nCorrect usage: hwe infile outfile.\n");
    exit_value = 1;
  }
  if ( ( *infile = fopen (argv[1], "r")) == (FILE *) NULL ) {
    REprintf ( "Can't read %s\n", argv[1]);
    exit_value = 2;
  }
  if( ( *outfile = fopen (argv[2], "w")) == (FILE *) NULL ) {
    REprintf ("Can't write %s\n", argv[2]);
    exit_value = 3;
  }
  return (exit_value);
}

int read_data (int *a, int *no_allele, int *total, randomization *sample, FILE **infile )
{
  register int i, j, l, err = 1;

  *total = 0;
  if( fscanf(*infile, "%d", no_allele) != 1) {
    REprintf("Please supply number of alleles\n");
    return ( err );
  }
  if ( *no_allele < 3 ) {
    REprintf("***Error! Number of alleles less than 3. \n");
    return ( err );
  }
  for ( i = 0; i < *no_allele; ++i ) {
    for ( j = 0; j <= i; ++j ) {
         l = LL(i, j);
         fscanf (*infile, "%d ", &a[l]);
         *total += a[l];
    }
  }
  if( fscanf(*infile, "%d %d %d \n", &sample->step,
                   &sample->group, &sample->size) != 3 ) {
    REprintf( " Please supply parameters.\n" );
    return ( err );
  }
  if ( sample->step < 1 || sample->group <= 1 ) {
    REprintf( "***Error in parameter specification.\n" );
    return ( err );
  }
  return (0);
}

void print_data (int *a, int no_allele, randomization sample, FILE **outfile)
{
  register int i, j, k, l;
  char line[250];

  line[0] = '-';
  k = 1;
  fprintf (*outfile, "Observed genotype frequencies: \n\n");
  for ( i = 0; i < no_allele; ++i ) {
    for ( j = k; j < k + 5; ++j ) line[j] = '-';
    line[j] = STR_END;
    k = j;
    fprintf (*outfile, "%s\n", line);
    fprintf (*outfile, "|");
    for ( j = 0; j <= i; ++j ) {
         l = LL(i, j);
         fprintf(*outfile, "%4d|", a[l]);
    }
    fprintf (*outfile, "\n");
  }
  fprintf (*outfile, "%s\n\n", line);
  fprintf (*outfile, "Total number of alleles: %2d\n\n", no_allele);
  fprintf(*outfile, "Number of initial steps: %d\n", sample.step);
  fprintf(*outfile, "Number of chunks: %d\n", sample.group);
  fprintf(*outfile, "Size of each chunk: %d\n\n", sample.size);
}

double cal_prob (int *a, Index index, double ln_p_old, int *actual_switch)
{
  double p1_ratio, p2_ratio;
  register double ln_p_new;
  double rand_num;
  int switch_ind, type;

  *actual_switch = 0;
  test_switch( a, index, &switch_ind, &type, &p1_ratio, &p2_ratio);
  switch (switch_ind)
    {
    case 0:
      ln_p_new = ln_p_old;
      break;
    case 1:
      if ( type == 1 ) p1_ratio = p2_ratio;
      rand_num = unif_rand();
      if ( rand_num < TRANS( p1_ratio ) ) {
           ndo_switch ( a, index, type );
           ln_p_new = ln_p_old + log (p1_ratio);
           *actual_switch = 1;
      } else ln_p_new = ln_p_old;
      break;
    default:
      rand_num = unif_rand();
      if ( rand_num <= TRANS(p1_ratio)) {
         ndo_switch ( a, index, 0 );
         ln_p_new = ln_p_old + log (p1_ratio);
         *actual_switch = 2;
       } else if ( rand_num <= TRANS(p1_ratio) + TRANS(p2_ratio) ) {
         ndo_switch ( a, index, 1 );
         ln_p_new = ln_p_old + log (p2_ratio);
         *actual_switch = 2;
       } else ln_p_new = ln_p_old;
       break;
    }
  return (ln_p_new);
}

void test_switch (int *a, Index index, int *switch_ind, int *switch_type, double *p1_rt, double *p2_rt)
{
  int k11, k22, k12, k21;

  *switch_ind = 0;
  k11 = L(index.i1, index.j1);
  k22 = L(index.i2, index.j2);
  k12 = L(index.i1, index.j2);
  k21 = L(index.i2, index.j1);
  if ( index.type <= 1 ) {
    if ( a[k11] > 0 && a[k22] > 0 ) {
         *switch_ind = 1;
         *switch_type = 0;
         *p1_rt = RATIO(a[k11], a[k12]) *  RATIO(a[k22], a[k21]) * index.cst;
    }
    if ( a[k12] > 0 && a[k21] > 0 ) {
         *switch_ind += 1;
         *switch_type = 1;
      *p2_rt = RATIO(a[k12], a[k11]) *  RATIO(a[k21], a[k22]) / index.cst;
    }
  } else {
    if ( a[k11] > 0 && a[k22] > 0 ) {
      *switch_ind = 1;
      *switch_type = 0;
         *p1_rt = RATIO(a[k11],a[k12] + 1.0)*RATIO(a[k22],a[k12]) * index.cst;
    }
    if ( a[k12] > 1 ) {
      *switch_ind += 1;
      *switch_type = 1;
      *p2_rt = RATIO(a[k12],a[k11]) * RATIO(a[k12] - 1,a[k22]) / index.cst;
    }
  }
}

#endif
