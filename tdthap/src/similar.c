#include <math.h>
#include <stdlib.h>

int nloci_sav;
double *spacing_sav = (double *) 0;
int focus_sav;
double power_sav;

/* S/R callable */

void set_tdt_similarity(int *nloci, double *spacing, int *focus, double *power) {
  int i;
  if (spacing_sav) free(spacing_sav);
  nloci_sav = *nloci;
  spacing_sav = (double *) calloc(nloci_sav+1, sizeof(double));
  for (i=0; i<=nloci_sav; i++) spacing_sav[i] = spacing[i];
  focus_sav = *focus;
  power_sav = *power;
}

void get_tdt_similarity(int *nloci, double *spacing, int *focus, double *power) {
  int i;
  *nloci = nloci_sav;
  for (i=0; i<=nloci_sav; i++) spacing[i] = spacing_sav[i];
  *focus = focus_sav;
  *power = power_sav;
}

double tdt_similarity(int *a, int *b) {
  double s;
  int i;
  /* Haplotypes must have focal locus in common */
  if (a[focus_sav-1] != b[focus_sav-1]) return(0.0);
  /* Find beginning of common area */
  for (i=focus_sav-2; i>=0 && a[i]==b[i]; i--) {}
  i++;
  s = spacing_sav[i++]/2.0; /* off-end correction */
  /* Find end of common area */
  for (; i<nloci_sav && a[i]==b[i]; i++) s += spacing_sav[i];
  s += spacing_sav[i]/2.0; /* off-end correction */
  if (power_sav!=1.0) s = pow(s, power_sav);
  return(s);
}
