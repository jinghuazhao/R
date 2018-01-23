void score_pairs(int *allele, int *nn, double *nplscore)
{
  int i, j, n = *nn;
  double score, total;

  *nplscore = 0.0;
  if (n < 2) return;
  score = total = 0.0;

  for (i=0; i<n; ++i)
    for (j=i+1; j<n; ++j) {
        if (allele[i*2] == allele[j*2]) score += 0.25;
        if (allele[i*2] == allele[j*2+1]) score += 0.25;
        if (allele[i*2+1] == allele[j*2]) score += 0.25;
        if (allele[i*2+1] == allele[j*2+1]) score += 0.25;
        total += 1.0;
    }

  *nplscore = (score/total);
}

void score_all(int *allele, int *nn, double *nplscore)
{
  int choices, i, j, a[3], bit[3], hn(int*,int);
  int n = *nn;
  double score;

  choices = 1;
  for (i=0; i<n; ++i)
    choices *= 2;

  score = 0;
  for (i=0; i<choices; ++i)
  {
    for (j=0; j<n; ++j) bit[j] = (i >> j) % 2;
    for (j=0; j<n; ++j) a[j] = allele[j*2+bit[j]];

    score += hn(a,n);
  }
  *nplscore = (score/(double)choices);
}

int hn(int *all, int n)
{
  int i, k, cnt, flg, seen[3], num[3], perms;
  cnt = 0;
  for (i=0; i<n; ++i)
  {
    k=0;
    flg = 1;
    while (k < cnt)
    {
      if (all[i] == seen[k])
      {
        flg =0;
        ++num[k];
      }
      ++k;
    }
    if (flg)
    {
      seen[cnt] = all[i];
      num[cnt] = 1;
      ++cnt;
    }
  }
  perms = 1;
  for (i=0; i<cnt; ++i)
    for (k=1; k<=num[i]; ++k)
      perms *= k;

  return perms-1;
}
