#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int    Y1[100];         /* frequency of case */
int    Y2[100];         /* frequency of control */
int    Y[100];          /* Y1 + Y2 */
int    N1;              /* sum Y1[] */
int    N2;              /* sum Y2[] */
int    N;               /* N1 + N2 */
int    K;               /* num of category */
int    Z[100];          /* sum_1-j Y1[i] */
int    S[100];          /* sum_1-j Y[i] */
double C_obs;           /* observed max chi^2 (among all 1 to the others) */
int    maxcol_obs;      /* column number of observed max chi^2  */
char   line[512][1000];
double Cout_obs;        /* observed max chi^2 (among all 1 to the others) */
double C_max_obs;       /* max(C_obs, Cout_obs) */
int    Ccol_obs;        /* column number for observed max */
int    Coutcol_obs;     /* column number for observed max */
int    Chi2Flag;        /* max chi2 is larger    --> 1 */
                        /* 1 to others is larger --> 0 */
double Combi(int, int);
void BasicStatistic(), CheckZero();
int CalcLj(int), CalcUj(int), CalcLoutj(int),CalcUoutj(int);
double MaxChiSquare(), MaxAmongOneToOthers();

void x22k(int *a, int *tablen, double *x2a, double *x2b, int *col1, int *col2, double *p)
{
    int      i, j;
    int      Lz[100], Uz[100];
    int      Lzsum,Uzsum;
    double   F[5000], F_prev[5000];
    int      Loop;
    int      TmpInt;
    int      ColumnONE = 100;
    int      ColumnACCUM = 100;
    int      CutFlag = 0;

    if(ColumnONE < 100) CutFlag = 1;
    if(ColumnACCUM < 100) CutFlag = 2;

    K=*tablen;
    for(i=0;i<K;i++)
    {
      Y1[i]=a[i];
      Y2[i]=a[i+K];
    }
    BasicStatistic();
    CheckZero();
    Rprintf("\n");
    Rprintf("Data : \n         ");
    for(i = 0; i < K; i++) Rprintf("%3d ",i + 1);
    Rprintf("Total\n-----------------------------------\nCase     ");
    for(i = 0; i < K; i++) Rprintf("%3d ",Y1[i]);
    Rprintf("%4d\nControl  ",N1);
    for(i = 0; i < K; i++) Rprintf("%3d ",Y2[i]);
    Rprintf("%4d\n-----------------------------------\nTotal    ",N2);
    for(i = 0; i < K; i++) Rprintf("%3d ",Y[i]);
    Rprintf("%4d\n",N);
    Rprintf("-------------------------------------------------------\n");

    C_obs = MaxChiSquare();
    Cout_obs = MaxAmongOneToOthers();
    *x2a=Cout_obs;
    *col1=Ccol_obs+1;
    *x2b=C_obs;
    *col2=Coutcol_obs+1;
    if(C_obs > Cout_obs) {
        C_max_obs = C_obs;
        maxcol_obs = Ccol_obs;
        Chi2Flag = 1;
        Rprintf("observed max is %f\n",C_max_obs);
        Rprintf("max chi2 %d...%d vs %d...%d\n",1,maxcol_obs+1,maxcol_obs+2,K);
        Rprintf("-------------------------------------------\n");
    }
    else {
        C_max_obs = Cout_obs;
        maxcol_obs = Coutcol_obs;
        Chi2Flag = 0;
        Rprintf("observed max is %f\n",C_max_obs);
        Rprintf("max 1-to-thers %d vs off %d\n",maxcol_obs + 1,maxcol_obs + 1);
        Rprintf("-------------------------------------------\n");
    }
    if(CutFlag == 1) {
      maxcol_obs = ColumnONE - 1;
      Chi2Flag = 0;
      Rprintf("calculate Pr(T >= max(%d vs others))\n",maxcol_obs + 1);
    }
    if(CutFlag == 2) {
      maxcol_obs = ColumnACCUM-1;
      Chi2Flag = 1;
      Rprintf("calculate Pr(T >= max(-%d vs %d-))\n",maxcol_obs+1,maxcol_obs+2);
    }

    Lz[0] = CalcLj(0);
    TmpInt = CalcLoutj(0);
    if(Lz[0] < TmpInt) Lz[0] = TmpInt;
    if(Lz[0] < 0) Lz[0] = 0;
    if(Lz[0] < Y[0] - N2) Lz[0] = Y[0] - N2;
    Uz[0] = CalcUj(0);
    TmpInt = CalcUoutj(0);
    if(Uz[0] > TmpInt) Uz[0] = TmpInt;
    if(Uz[0] > Y[0]) Uz[0] = Y[0];
    if(Uz[0] > N1) Uz[0] = N1;
    for(i = 1; i < K - 1; i++)
    {
        Lz[i] = CalcLj(i);
        TmpInt = CalcLoutj(i) + Lz[i-1];
        if(Lz[i] < TmpInt) Lz[i] = TmpInt;
        if(Lz[i] < 0) Lz[i] = 0;
        if(Lz[i] < S[i] - N2) Lz[i] = S[i] - N2;
        Uz[i] = CalcUj(i);
        TmpInt = CalcUoutj(i) + Uz[i-1];
        if(Uz[i] > TmpInt) Uz[i] = TmpInt;
        if(Uz[i] > N1) Uz[i] = N1;
        if(Uz[i] > S[i]) Uz[i] = S[i];
    }
    for(i = Lz[0]; i <= Uz[0]; i++) F_prev[i] = 1.0;
    Loop = 1;
    while(1)
    {
        for(i = Lz[Loop]; i <= Uz[Loop]; i++)
        {
            F[i] = 0.0;
            Lzsum = CalcLj(Loop - 1);
            if(Lzsum < 0) Lzsum = 0;
            if(Lzsum < i - Y[Loop]) Lzsum = i - Y[Loop];
            TmpInt = i - CalcUoutj(Loop);
            if(Lzsum < TmpInt) Lzsum = TmpInt;
            Uzsum = CalcUj(Loop - 1);
            if(Uzsum > S[Loop]) Uzsum = S[Loop];
            if(Uzsum > i) Uzsum = i;
            TmpInt = i - CalcLoutj(Loop);
            if(Uzsum > TmpInt) Uzsum = TmpInt;
            for(j = Lzsum; j <= Uzsum; j++)  F[i] += F_prev[j]
                  * Combi(i, j) * Combi(S[Loop] - i, S[Loop-1] - j)
                    / Combi(S[Loop], S[Loop-1]);
            if(F[i] > 1.0) F[i] = 1.0;
            if(F[i] < 0.0) F[i] = 0.0;
        }
        for(i = Lz[Loop - 1]; i <= Uz[Loop - 1]; i++) F_prev[i] = 0.0;
        for(i = Lz[Loop]; i <= Uz[Loop]; i++){
            F_prev[i] = F[i];
            F[i] = 0.0;
        }
        if(Loop == (K - 2)) break;
        Loop++;
    }
    F[N1] = 0.0;
    Lzsum = CalcLj(K-2);
    if(Lzsum < 0) Lzsum = 0;
    if(Lzsum < N1 - Y[K-1]) Lzsum = N1 - Y[K-1];
    TmpInt = N1 - CalcUoutj(K-1);
    if(Lzsum < TmpInt) Lzsum = TmpInt;
    Uzsum = CalcUj(K-2);
    if(Uzsum > S[K-2]) Uzsum = S[K-2];
    if(Uzsum > N1) Uzsum = N1;
    TmpInt = N1 - CalcLoutj(K-1);
    if(Lzsum > TmpInt) Lzsum = TmpInt;
    for(i = Lzsum; i <= Uzsum; i++) F[N1] += F_prev[i]
          * Combi(N1, i) * Combi(N2, S[K - 2] - i) / Combi(N, S[K - 2]);
    *p = 1.0 - F[N1];
    Rprintf("p-value = %.10f\n",1.0 - F[N1]);
    Rprintf("------------------------------------------------\n");
}

double MaxChiSquare()
{
  int     i,j;
  int     z11,z12,z21,z22;
  int     c1,c2;
  double  chi[100];
  double  maxchi;
  int     maxindex=0;
  int     z11WM=0,z12WM=0,z21WM=0,z22WM=0;
  int     c1WM=0,c2WM=0;
  maxchi = -1.0;
  Rprintf("# ------------------------------------\n");
  Rprintf("# cut point   statistic values\n");
  for(i = 0; i < K - 1; i++)
  {
    z11 = 0;
    z12 = 0;
    z21 = 0;
    z22 = 0;
    c1 = 0;
    c2 = 0;
    for(j = 0; j <= i; j++)
    {
      z11 += Y1[j];
      z21 += Y2[j];
      c1 += Y[j];
    }
    for(j = i + 1; j < K; j++)
    {
      z12 += Y1[j];
      z22 += Y2[j];
      c2 += Y[j];
    }
    chi[i] = (double)N *
      (double)((z11 * z22) - (z12 * z21)) * ((z11 * z22) - (z12 * z21))
      / ((double)N1 * N2 * c1 * c2);
    Rprintf("#    %d-%d      %f\n",i+1,i+2,chi[i]);
    if(maxchi < chi[i]){
      maxchi = chi[i];
      maxindex = i;
      z11WM = z11;
      z12WM = z12;
      z21WM = z21;
      z22WM = z22;
      c1WM = c1;
      c2WM = c2;
    }
  }
  Rprintf("# ------------------------------------\n");
  Ccol_obs = maxindex;
  Rprintf("Max chi square = %f\n",maxchi);
  Rprintf("where the table is divided between\n");
  Rprintf("before the %d th and after the %d th category\n\n",maxindex + 1, maxindex + 2);
  Rprintf("1,...,%d  %d,...,%d\n",maxindex + 1, maxindex + 2,K);
  Rprintf("  %3d      %3d       %3d\n",z11WM,z12WM,N1);
  Rprintf("  %3d      %3d       %3d\n",z21WM,z22WM,N2);
  Rprintf("  %3d      %3d       %3d\n",c1WM,c2WM,N);
  Rprintf("-------------------------------------------------------\n");
  return maxchi;
}

double MaxAmongOneToOthers()
{
  int     i;
  int     z11,z12,z21,z22;
  int     c1,c2;
  double  maxchi;
  double  chi;
  int     maxindex=0;
  int     z11WM=0,z12WM=0,z21WM=0,z22WM=0;
  int     c1WM=0,c2WM=0;
  double  E1,E2;
  maxchi = -1.0;
  Rprintf("# ----------------------------\n");
  Rprintf("# considered column   statistic values\n");
  for(i = 0; i < K; i++)
  {
    z11 = Y1[i];
    z12 = N1 - Y1[i];
    z21 = Y2[i];
    z22 = N2 - Y2[i];
    c1 = Y[i];
    c2 = N - Y[i];
    E1 = (double) c1 * N1 / N;
    E2 = (double) c1 * N2 / N;

    chi = (double)N *
      (double)((z11 * z22) - (z12 * z21)) * ((z11 * z22) - (z12 * z21))
      / ((double)N1 * N2 * c1 * c2);
    Rprintf("#    %d                %f\n",i+1,chi);
    if(chi > maxchi)
    {
      maxchi = chi;
      maxindex = i;
      z11WM = z11;
      z12WM = z12;
      z21WM = z21;
      z22WM = z22;
      c1WM = c1;
      c2WM = c2;
    }
  }
  Rprintf("# ----------------------------\n");
  Coutcol_obs = maxindex;
  Rprintf("Max Chi Square (among all 1-to-others) = %f\n",maxchi);
  Rprintf("where the table is divided between\n");
  Rprintf("the %d th category and the others\n\n",maxindex + 1);
  Rprintf(" %d th  the others\n",maxindex + 1);
  Rprintf("  %3d      %3d       %3d\n",z11WM,z12WM,N1);
  Rprintf("  %3d      %3d       %3d\n",z21WM,z22WM,N2);
  Rprintf("  %3d      %3d       %3d\n",c1WM,c2WM,N);
  Rprintf("-------------------------------------------------------\n");
  return maxchi;
}

double Combi(int x, int y)
{
  int i;
  double z;
  if(y == 0) z = 1.0;
  else {
    if((double)y > 0.5 * x) y = x - y;
    z = 1.0;
    for(i = 0; i < y; i++){
      z *= ((x - (double)i) / (y - (double)i));
    }
  }
  return z;
}

void BasicStatistic()
{
  int i;
  N1 = 0;
  N2 = 0;
  N = 0;
  for(i = 0; i < K; i++) 
  {
    Y[i] = Y1[i] + Y2[i];
    N1 += Y1[i];
    N2 += Y2[i];
  }
  N = N1 + N2;
}

void CheckZero()
{
  int i,j;
  for(i = 0; i < K; i++)
  {
    if(Y[i] == 0)
    {
      K--;
      for(j = i; j < K; j++)
      {
        Y1[j] = Y1[j + 1];
        Y2[j] = Y2[j + 1];
        Y[j] = Y[j + 1];
      }
      i--;
    }
    Z[i] = 0;
    S[i] = 0;
  }
  Z[0] = Y1[0];
  S[0] = Y[0];
  for(i = 1; i < K; i++)
  {
    Z[i] = Z[i - 1] + Y1[i];
    S[i] = S[i - 1] + Y[i];
  }
}

int CalcLj(int col)
{
  double L;
  double a,b;
  int    Lint;
  a = (double)N1 * S[col] / N;
  if(Chi2Flag) 
  {
    /* max chi2 is larger */
    b = (double)(1.0 / N) * sqrt((double)S[col] * (N - S[col]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) / (S[maxcol_obs] * (N - S[maxcol_obs])));
  }
  else 
  {
    /* outlier type is larger */
    b = (double)(1.0 / N) * sqrt((double)S[col] * (N - S[col]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) / (Y[maxcol_obs] * (N - Y[maxcol_obs])));
  }
  L = a - b;
  Lint = (int)floor(L + 1.0);
  return Lint;
}

int CalcUj(int col)
{
  double U;
  double a,b;
  int    Uint;
  a = (double)N1 * S[col] / N;
  if(Chi2Flag)
  {
    /* max chi2 is larger */
    b = (double)(1.0 / N) * sqrt((double)S[col] * (N - S[col]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) / (S[maxcol_obs] * (N - S[maxcol_obs])));
  }
  else
  {
    /* outlier type is larger */
    b = (double)(1.0 / N) * sqrt((double)S[col] * (N - S[col]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) / (Y[maxcol_obs] * (N - Y[maxcol_obs])));
  }
  U = a + b;
  Uint = (int)ceil(U - 1.0);
  return Uint;
}

int CalcLoutj(int col)
{
  double L;
  double a,b;
  int    Lint;

  a = (double)N1 * Y[col] / N;
  if(Chi2Flag)
  {
    /* max chi2 is larger */
    b = (double)(1.0 / N) * sqrt((double)Y[col] * (N - Y[col]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) / (S[maxcol_obs] * (N - S[maxcol_obs])));
  }
  else
  {
    /* outlier type is larger */
    b = (double)(1.0 / N) * sqrt((double)Y[col] * (N - Y[col]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) / (Y[maxcol_obs] * (N - Y[maxcol_obs])));
  }
  L = a - b;
  Lint = (int)floor(L + 1.0);
  return Lint;
}

int CalcUoutj(int col)
{
  double U;
  double a,b;
  int    Uint;

  a = (double)N1 * Y[col] / N;
  if(Chi2Flag)
  {
    /* max chi2 is larger */
    b = (double)(1.0 / N) * sqrt((double)Y[col] * (N - Y[col]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) * (N * Z[maxcol_obs] - N1 * S[maxcol_obs]) / (S[maxcol_obs] * (N - S[maxcol_obs])));
  }
  else
  {
    /* outlier type is larger */
    b = (double)(1.0 / N) * sqrt((double)Y[col] * (N - Y[col]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) * (N * Y1[maxcol_obs] - N1 * Y[maxcol_obs]) / (Y[maxcol_obs] * (N - Y[maxcol_obs])));
  }
  U = a + b;
  Uint = (int)ceil(U - 1.0);
  return Uint;
}

#ifdef executable
void ReadData()
{
  int i;
    while(1){
      Rprintf("num of category : ");
      fgets(line[0],sizeof(line[0]),stdin);
      sscanf(line[0],"%d",&K);
      if(K > 100) Rprintf("Sorry, num of category must be less than 100\n");
      else break;
    }
    for(i = 0; i < K; i++){
      Rprintf("N[1,%d] : ",i+1);
      fgets(line[i],sizeof(line[i]),stdin);
      sscanf(line[i],"%d",&Y1[i]);
    }
    for(i = 0; i < K; i++){
      Rprintf("N[2,%d] : ",i+1);
      fgets(line[i],sizeof(line[i]),stdin);
      sscanf(line[i],"%d",&Y2[i]);
    }
}

void DisplayData()
{
  int i;
  Rprintf("\n");
  Rprintf("Data : \n         ");
  for(i = 0; i < K; i++){
    Rprintf("%3d ",i + 1);
  }
  Rprintf("Total\n-----------------------------------\nCase     ");
  for(i = 0; i < K; i++){
    Rprintf("%3d ",Y1[i]);
  }
  Rprintf("%4d\nControl  ",N1);
  for(i = 0; i < K; i++){
    Rprintf("%3d ",Y2[i]);
  }
  Rprintf("%4d\n-----------------------------------\nTotal    ",N2);
  for(i = 0; i < K; i++){
    Rprintf("%3d ",Y[i]);
  }
  Rprintf("%4d\n",N);
  Rprintf("-------------------------------------------------------\n");
}

int main(int argc,char **argv)
{
    int      i, j;
    int      Lz[100], Uz[100];
    int      Lzsum,Uzsum;
    double   F[5000], F_prev[5000];
    int      Loop;
    int      TmpInt;
    int      ColumnONE = 100;
    int      ColumnACCUM = 100;
    int      CutFlag = 0;

    while((argc > 1) && (argv[1][0] == '-')){
      switch(argv[1][1]){
      case 'C':
        ColumnONE = atoi(&argv[1][2]); break;
      case 'A':
        ColumnACCUM = atoi(&argv[1][2]); break;
      }
      argv++;
      argc--;
    }
    if(ColumnONE < 100) CutFlag = 1;
    if(ColumnACCUM < 100) CutFlag = 2;

    ReadData();
    BasicStatistic();
    CheckZero();
    DisplayData();

    C_obs = MaxChiSquare();
    Cout_obs = MaxAmongOneToOthers();
    if(C_obs > Cout_obs){
        C_max_obs = C_obs;
        maxcol_obs = Ccol_obs;
        Chi2Flag = 1;
        Rprintf("observed max is %f\n",C_max_obs);
        Rprintf("max chi2 %d...%d vs %d...%d\n",1,maxcol_obs+1,maxcol_obs+2,K);
        Rprintf("-------------------------------------------\n");
    }
    else{
        C_max_obs = Cout_obs;
        maxcol_obs = Coutcol_obs;
        Chi2Flag = 0;
        Rprintf("observed max is %f\n",C_max_obs);
        Rprintf("max 1-to-thers %d vs off %d\n",maxcol_obs + 1,maxcol_obs + 1);
        Rprintf("-------------------------------------------\n");
    }
    if(CutFlag == 1){
      maxcol_obs = ColumnONE - 1;
      Chi2Flag = 0;
      Rprintf("calculate Pr(T >= max(%d vs others))\n",maxcol_obs + 1);
    }
    if(CutFlag == 2){
      maxcol_obs = ColumnACCUM-1;
      Chi2Flag = 1;
      Rprintf("calculate Pr(T >= max(-%d vs %d-))\n",maxcol_obs+1,maxcol_obs+2);
    }

    Lz[0] = CalcLj(0);
    TmpInt = CalcLoutj(0);
    if(Lz[0] < TmpInt) Lz[0] = TmpInt;
    if(Lz[0] < 0) Lz[0] = 0;
    if(Lz[0] < Y[0] - N2) Lz[0] = Y[0] - N2;
    Uz[0] = CalcUj(0);
    TmpInt = CalcUoutj(0);
    if(Uz[0] > TmpInt) Uz[0] = TmpInt;
    if(Uz[0] > Y[0]) Uz[0] = Y[0];
    if(Uz[0] > N1) Uz[0] = N1;
    for(i = 1; i < K - 1; i++){
        Lz[i] = CalcLj(i);
        TmpInt = CalcLoutj(i) + Lz[i-1];
        if(Lz[i] < TmpInt) Lz[i] = TmpInt;
        if(Lz[i] < 0) Lz[i] = 0;
        if(Lz[i] < S[i] - N2) Lz[i] = S[i] - N2;
        Uz[i] = CalcUj(i);
        TmpInt = CalcUoutj(i) + Uz[i-1];
        if(Uz[i] > TmpInt) Uz[i] = TmpInt;
        if(Uz[i] > N1) Uz[i] = N1;
        if(Uz[i] > S[i]) Uz[i] = S[i];
    }
    for(i = Lz[0]; i <= Uz[0]; i++){
        F_prev[i] = 1.0;
    }
    Loop = 1;
    while(1){
        for(i = Lz[Loop]; i <= Uz[Loop]; i++){
            F[i] = 0.0;
            Lzsum = CalcLj(Loop - 1);
            if(Lzsum < 0) Lzsum = 0;
            if(Lzsum < i - Y[Loop]) Lzsum = i - Y[Loop];
            TmpInt = i - CalcUoutj(Loop);
            if(Lzsum < TmpInt) Lzsum = TmpInt;
            Uzsum = CalcUj(Loop - 1);
            if(Uzsum > S[Loop]) Uzsum = S[Loop];
            if(Uzsum > i) Uzsum = i;
            TmpInt = i - CalcLoutj(Loop);
            if(Uzsum > TmpInt) Uzsum = TmpInt;
            for(j = Lzsum; j <= Uzsum; j++){
                F[i] += F_prev[j]
                  * Combi(i, j) * Combi(S[Loop] - i, S[Loop-1] - j)
                    / Combi(S[Loop], S[Loop-1]);
            }
            if(F[i] > 1.0) F[i] = 1.0;
            if(F[i] < 0.0) F[i] = 0.0;
        }
        for(i = Lz[Loop - 1]; i <= Uz[Loop - 1]; i++){
            F_prev[i] = 0.0;
        }
        for(i = Lz[Loop]; i <= Uz[Loop]; i++){
            F_prev[i] = F[i];
            F[i] = 0.0;
        }
        if(Loop == (K - 2)) break;
        Loop++;
    }
    F[N1] = 0.0;
    Lzsum = CalcLj(K-2);
    if(Lzsum < 0) Lzsum = 0;
    if(Lzsum < N1 - Y[K-1]) Lzsum = N1 - Y[K-1];
    TmpInt = N1 - CalcUoutj(K-1);
    if(Lzsum < TmpInt) Lzsum = TmpInt;
    Uzsum = CalcUj(K-2);
    if(Uzsum > S[K-2]) Uzsum = S[K-2];
    if(Uzsum > N1) Uzsum = N1;
    TmpInt = N1 - CalcLoutj(K-1);
    if(Lzsum > TmpInt) Lzsum = TmpInt;
    for(i = Lzsum; i <= Uzsum; i++){
        F[N1] += F_prev[i]
          * Combi(N1, i) * Combi(N2, S[K - 2] - i)
            / Combi(N, S[K - 2]);
    }
    Rprintf("p-value = %.10f\n",1.0 - F[N1]);
    Rprintf("------------------------------------------------\n");
}

/******************************************************************
 *     program for calculating p values of maximum of
 *     max all one-to-others and max accumulated
 *     statistic for 2 x K contingency tables
 *
 *     option:
 *        -Ccolumn : calculate upper probability
 *                   for [column] vs [others]
 *        -Acolumn : calculate upper probability
 *                   for [1,..,column] vs [column+1,..,K]
 *
 *     last update : July, 12. 2002
*******************************************************************/
#endif

