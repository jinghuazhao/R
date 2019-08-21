#include <R.h>
#include <stdlib.h>
#include <stdio.h>

#include "nuc_fam.h"

/* 
   Declarations of functions local to this file
*/

int inherited(int allele, int p_gtype[2]); /* Possible inheritance of allele */
int poss_tr(int a_f, int a_m, int o_g[2]); /* Possible transmission pattern */
int trans(int phase[], int inhvec[], int m);  /* Phase and transmission */
int fill_in(int child[2], int unknown[2], int known[2]);

/* ============================ S/R-callable ================================*/

/* 
   Take a pedfile-style input and write to file the affected offspring 
   plus the transmitted and untransmitted haplotypes. First argument is 
   length of input vectors, but is returned as the number of records 
   written. Treatment of multiple cases per family is either:
   use all (0),
   use all but generate separate family for each case (with all other 
   sibs non-cases) (1), or
   use first case only (2).
   Imputation of missing parental genotypes may be done either with or 
   without use of affected offspring.
*/

void hap_transmit(int *n, int *ped, int *id, int *father, int *mother,
		  int *sex, int *aff, int *if_qt, double *qt, 
		  int *m, int *markers, 
		  int *multiple_cases, int *impute_using_affected,
		  char **ofname) {
  Family *first, *f, *prev;
  FILE *outfile;
  int nn, mm, hr, iqt;
  char *tmp;
  nn = *n;
  mm = *m;
  iqt = *if_qt;
  if (!*if_qt) qt = (double *) 0;
  first = nuclear(nn, ped, id, father, mother, sex, aff, qt, mm, markers);
  /* Multiple case treatment */
  if (*multiple_cases) {
    for (f=first; f; f=f->next) {
      if (*multiple_cases == 1) {
	prev = f;
	f = expand_family(f, mm);
	if (!f) goto overflow;
      }
      else if (*multiple_cases == 2) {
	use_only_first(f);
      }
    }
  }

  /* Do remaining computations on families */
  
  prev = (Family *) 0;
  for (f=first; f; f=f->next) {
    /* Impute missing parental genotypes */
    impute_parent(f, mm, (int) *impute_using_affected);
    /* Compute inheritance vectors */
    inheritance(f, mm);
    /* Resolve haplotype phase and transmission */
    hr = haplotype(f, mm, 1);
    /* If recombination, write error message */
    if (hr<0) {
      REprintf("*** Recombination/expaternity at locus %d *** ", -hr);
      show_family(f);
    }
    /* If no information or recombination, omit family */
    if (hr!=0) {
      if (prev)
	prev->next = f->next;
      else
	first = f->next;
    }
    else {
      prev = f;
    }
  } 

  /* Write haplotypes to disk */

  tmp = *ofname;
  /* If no file name supplied, generate one */
  if (!*tmp) {
    tmp = mkstemp((char *) 0);
    *ofname = tmp;
  }
  outfile = fopen(tmp, "wb");
  if (outfile) {
    *n = hap_write(first, mm, iqt, outfile);
    fclose(outfile);
  }
  else {
    REprintf("*** Couldn't open temporary file %s\n", tmp);
    *n = 0;
  }

  /* Now return memory to system */

  while (first) {
    f = first;
    first = first->next;
    del_family(f);
  }
  return;

  /* Memory overflow */

overflow:
  warn("Memory overflow while or after expanding family", f);

}

/*
  Read haplotypes back into int arrays, and delete file
*/

void hap_read(int *n, int *ped, int *id, int *father, int *mother,
	      int *if_qt, double *qt, 
	      int *m, int *f_tr, int *f_un, int *m_tr, int *m_un, 
	      char **ifname) {
  FILE *infile;
  int locus, mm, nn, i, ij, fr, n_read;
  int *hap;
  infile = fopen(*ifname, "rb");
  if (!infile) {
    *n = 0;
    return;
  }
  if (!*if_qt) qt = (double *) 0;
  nn = *n;
  mm = *m;
  hap = (int *) calloc(mm, sizeof(int));
  for (i=0; i<nn;  i++) {
    n_read = i;
    fr = fread(ped+i, sizeof(int), 1, infile);
    if (!fr) goto end_of_file;
    fr = fread(id+i, sizeof(int), 1, infile);
    if (!fr) goto end_of_file;
    fr = fread(father+i, sizeof(int), 1, infile);
    if (!fr) goto end_of_file;
    fr = fread(mother+i, sizeof(int), 1, infile);
    if (!fr) goto end_of_file;
    if (qt) {
      fr = fread(qt+i, sizeof(double), 1, infile);
      if (!fr) goto end_of_file;
    }
    fr = fread(hap, sizeof(int), mm, infile);
    if (fr<mm) goto end_of_file;
    for (ij=i, locus=0; locus<mm; locus++, ij += nn) f_tr[ij] = hap[locus];
    fr = fread(hap, sizeof(int), mm, infile);
    if (fr<mm) goto end_of_file;
    for (ij=i, locus=0; locus<mm; locus++, ij += nn) f_un[ij] = hap[locus];
    fr = fread(hap, sizeof(int), mm, infile);
    if (fr<mm) goto end_of_file;
    for (ij=i, locus=0; locus<mm; locus++, ij += nn) m_tr[ij] = hap[locus];
    fr = fread(hap, sizeof(int), mm, infile);
    if (fr<mm) goto end_of_file;
    for (ij=i, locus=0; locus<mm; locus++, ij += nn) m_un[ij] = hap[locus];
  }
  /* Normal termination */
  remove(*ifname);
  free(hap);
  return;
end_of_file:
  remove(*ifname);
  free(hap);
  *n = n_read;
  return;
}

/* ============================ General ================================*/

/* 
   Break linkage-style pedfile data into a list of nuclear families
   Note that elements of mother and father arrays are set to zero when
   subjects are placed as children (one can't be a child in more than one 
   family). Arrays are (long) ints because the function might be called 
   originally from R or S. 
*/

Family *nuclear(int n, int *ped, int *mem, int *father, int *mother,
		int *sex, int *aff_status, double *qt, 
		int m, int *markers) {
  int person, j, pj, twom, placed, prev_placed;
  int p_ped, p_fat, p_mot, p_mem;
  Family *res, *f, *f2add, *f_last=NULL;
  Offspring *child, *prev;
  twom = 2*m;
  res = (Family *) 0;
  /* Repeatedly pass through persons until no-one else can be placed */
  prev_placed = -1;
  placed = 0;
  while (placed > prev_placed) {
    prev_placed = placed;
    for (person=0; person<n; person++) {
      p_ped = ped[person];
      p_mem = mem[person];
      p_fat = father[person];
      p_mot = mother[person];
      /* If person is an offspring ... */
      if (p_fat && p_mot) { 
	/* See if this person is a child of an existing family */
	for (f2add=(Family *) 0, f=res; f; f= f->next) {
	  if ((p_ped==f->pedigree) && (p_fat==f->father_id) && 
	      (p_mot==f->mother_id)) {
	    f2add = f;
	    break;
	  } 
	}
	/* If not, create new family record */
	if (!f2add) {
	  f2add = new_family(m);
	  if (!f2add) goto overflow;
	  f2add->pedigree = p_ped;
	  f2add->father_id = p_fat;
	  f2add->mother_id = p_mot;
	  if (res) 
	    f_last->next = f2add;
	  else 
	    res = f2add;
	  f_last = f2add;
	}
	/* Add child into family */
	prev = (Offspring *) 0;
	for (child=f2add->children; child; child=child->next) prev = child;
	child = new_child(m);
	if (!child) goto overflow;
	if (prev)
	  prev->next = child;
	else
	  f2add->children = child;
	child->id = p_mem;
	child->affected = aff_status[person];
	child->qt = qt? qt[person]: 0.0;
	child->sex = sex? sex[person]: 0;
	for (pj=person, j=0; j<twom; j++, pj+=n) 
	  child->markers[j] = markers[pj];
	father[person] = mother[person] = 0;
	placed++; /* Person placed as an offspring */
      }
      /* Whether or not an offspring, anyone can be a parent ... */
      for (f=res; f; f= f->next) {
	if (p_ped==f->pedigree) {
	  /* If family missing father, is this him? */
	  if (!(f->check % 2) &&  (p_mem==f->father_id)) {
	    if (sex && sex[person]!=1) {
	      warn("Father not male", f);
	    }
	    else {
	      f->check++;
	      for (pj=person, j=0; j<twom; j++, pj+=n) 
		f->father[j] = markers[pj];
	    }
	  }
	  /* If family missing mother, is this her? */
	  if (!(f->check / 2) && (p_mem==f->mother_id)) {
	    if (sex && sex[person]!=2) {
	      warn("Mother not female", f);
	    }
	    else {
	      f->check+=2;
	      for (pj=person, j=0; j<twom; j++, pj+=n) 
		f->mother[j] = markers[pj];
	    }
	  }
	}
      }
    }
  }
  return res;
overflow:
  warn("Memory overflow while or after storing family", f_last);
  return NULL;
}

/* 
   Write a list of transmitted and untransmitted htypes for affected children
   Return number of records written
*/

int hap_write(Family *first, int m, int if_qt, FILE *stream) {
  Family *f;
  Offspring *ch;
  int i, locus, ia, which, twom;
  int *haps;
  twom = 2*m;
  haps = (int *) calloc(twom, sizeof(int));
  for (i=0, f=first; f; f=f->next) {
    for (ch=f->children; ch; ch=ch->next) {
      /* Generate record if affected with either transmission known */
      if (ch->affected==2 && (ch->f_tr || ch->m_tr)) {
	i++;
	fwrite(&f->pedigree, sizeof(int), 1, stream); /* Pedigree */
	fwrite(&ch->id, sizeof(int), 1, stream); /*  Id */
	fwrite(&f->father_id, sizeof(int), 1, stream); /* Father */
	fwrite(&f->mother_id, sizeof(int), 1, stream); /* Mother */
	if (if_qt) fwrite(&ch->qt, sizeof(double), 1, stream); /* QT value */
	/* Paternal haplotypes , transmitted and untransmitted */
	if (ch->f_tr) {
	  for (ia=locus=0; locus<m; locus++, ia+=2) {
	    if (f->phase_f[locus]) {
	      which = (ch->f_tr + f->phase_f[locus]) % 2;
	      if (!which) { /* Order as stored */
		haps[locus] = f->father[ia];
		haps[locus+m] = f->father[ia+1];
	      }
	      else { /* Swap order */
		haps[locus] = f->father[ia+1];
		haps[locus+m] = f->father[ia];
	      }
	    }
	    else {
	      haps[locus] = haps[locus+m] = 0;
	    }
	  }
	}
	else {
	  for (ia=0; ia<twom; ia++) haps[ia] = 0;
	}
	fwrite(haps, sizeof(int), twom, stream); 
	/* Maternal haplotypes, transmitted and untransmitted */
	if (ch->m_tr) {
	  for (ia=locus=0; locus<m; locus++, ia+=2) {
	    if (f->phase_m[locus]) {
	      which = (ch->m_tr + f->phase_m[locus]) % 2;
	      if (!which) { /* Order as stored */
		haps[locus] = f->mother[ia];
		haps[locus+m] = f->mother[ia+1];
	      }
	      else { /* Swap order */
		haps[locus] = f->mother[ia+1];
		haps[locus+m] = f->mother[ia];
	      }
	    }
	    else {
	      haps[locus] = haps[locus+m] = 0;
	    }
	  }
	}
	else {
	  for (ia=0; ia<twom; ia++) haps[ia] = 0;
	}
	fwrite(haps, sizeof(int), twom, stream); 
      }
    }
  }
  free(haps);
  return i;
}

/* 
   Compute inheritance vectors for all children in the sibship 
   Returns n if error/ex-paternity for n-th sib, 0 otherwise
*/

int inheritance(Family *f, int m) {
  int ierror, locus, ia1, ia2,  
    f1, f2, m1, m2, i11, i12, i21, i22, ni;
  int *ch;
  Offspring *child;

  for (ierror=1, child= f->children; child; child = child->next, ierror++) {
    for (ia1=locus=0; locus<m; locus++, ia1+=2) {
      ia2 = ia1+1;
      f1 = f->father[ia1];
      f2 = f->father[ia2];
      m1 = f->mother[ia1];
      m2 = f->mother[ia2];
      ch = child->markers+ia1;
      /* Check possibility of four transmission patterns */
      i11 = poss_tr(f1, m1, ch);
      i12 = poss_tr(f1, m2, ch);
      i21 = poss_tr(f2, m1, ch);
      i22 = poss_tr(f2, m2, ch);
      ni = i11+i12+i21+i22;
      /* Carry on if any transmission patterns possible */
      if (ni) {
	/* If a unique pattern is selected, inheritance can be computed */
	if (ni==1) {
	  if (i11) {
	    child->ivec_f[locus] = 1;
	    child->ivec_m[locus] = 1;
	  }
	  else if (i12) {
	    child->ivec_f[locus] = 1;
	    child->ivec_m[locus] = 2;
	  }
	  if (i21) {
	    child->ivec_f[locus] = 2;
	    child->ivec_m[locus] = 1;
	  }
	  if (i22) {
	    child->ivec_f[locus] = 2;
	    child->ivec_m[locus] = 2;
	  }
	}
	/* else if there are two possibilities, one component might be */
	else if (ni==2) {
	  if (i11 && i12) 
	    child->ivec_f[locus] = 1;
	  else if (i21 && i22) 
	    child->ivec_f[locus] = 2;
	  else if (i11 && i21) 
	    child->ivec_m[locus] = 1;
	  else if (i12 && i22) 
	    child->ivec_m[locus] = 2;
	}
      }
    }
  }
  return 0;
}

/*
  Compute parental haplotypes and transmission
  Return 1 if uninformative, -L if recombination at locus L, 0 if successful
  If resolve_homozygous is True, we resolve the phase of the parental 
  haplotype arbitrarily (since it makes no difference), but transmission 
  remains unresolved.
*/

int haplotype(Family *f, int m, int resolve_homozygous) {
  int locus, nin_f, nin_m, max_f, max_m, ft, mt, ia1, ia2;
  int *ph_f, *ph_m, *iv_f, *iv_m;
  Offspring *child;

  /* Initialisation (just in case of uninformative families) */

  ph_f = f->phase_f;
  ph_m = f->phase_m;
  for (locus=0; locus<m; locus++) ph_f[locus] = ph_m[locus] = 0;

  /* On first pass, set phase vectors to most complete inheritance vectors */
  
  for (max_f=max_m=0, child = f->children; child; child = child->next) {
    iv_f = child->ivec_f;
    iv_m = child->ivec_m;
    for (nin_f=nin_m=locus=0;  locus<m; locus++) {
      if (iv_f[locus]) nin_f++;
      if (iv_m[locus]) nin_m++;
    }
    if (nin_f>max_f) {
      max_f = nin_f;
      for (locus=0; locus<m; locus++) ph_f[locus] = iv_f[locus];
    }
    if (nin_m>max_m) {
      max_m = nin_m;
      for (locus=0; locus<m; locus++) ph_m[locus] = iv_m[locus];
    }
  }
  if (!max_f && !max_m) return 1;

  /* 
     Repeatedly cycle through children until resolution not improved 
     max_f, max_m are number of loci resolved in the last pass
     nin_f, nin_m are the numbers resolved the time before
  */

  nin_f = nin_m = 0; /* to ensure loop gets executed first time */
  while (max_f>nin_f || max_m>nin_m) {
    for (child = f->children; child; child = child->next) {
      iv_f = child->ivec_f;
      iv_m = child->ivec_m;
      ft = trans(ph_f, iv_f, m);
      if (ft<0) return ft;  /* Recombination */
      mt = trans(ph_m, iv_m, m);
      if (mt<0) return mt;  /* Recombination */
      child->f_tr = ft;
      child->m_tr = mt;
    }
    nin_f = max_f;
    nin_m = max_m;
    /* How many resolved now? */
    for (max_f=locus=0; locus<m; locus++) if (ph_f[locus]) max_f++;
    for (max_m=locus=0; locus<m; locus++) if (ph_f[locus]) max_m++;
  }

  /* Optionally and arbitrarily resolve phase of homozygous loci */

  if (resolve_homozygous) {
    for (ia2=1, ia1=locus=0; locus<m; locus++, ia1+=2, ia2+=2) {
      if (f->father[ia1]==f->father[ia2]) {
	max_f++;
	ph_f[locus] = 1;
      }
      if (f->mother[ia1]==f->mother[ia2]) {
	max_m++;
	ph_m[locus] = 1;
      }
    }
  }

  return 0;
}

/* 
   Impute missing parental data when possible 
   return 1 if error/ex-paternity, 0 otherwise
   If use_affected is true, we use all offspring to fill in parental 
   alleles. Otherwise we use only unaffected and unknown offspring.
*/

int impute_parent(Family *f, int m, int use_affected) {
  int locus, ia1, ia2, miss_f, miss_m;
  Offspring *child;
  for (ia2=1, ia1=locus=0; locus<m; locus++, ia1+=2, ia2+=2) {
    miss_f = !f->father[ia1] || !f->father[ia2];
    miss_m = !f->mother[ia1] || !f->mother[ia2];
    /* If one parent only missing, we may be able to guess the rest */
    if ((miss_f && !miss_m)||(miss_m && !miss_f)) {
      for(child = f->children; child; child=child->next) {
	if (use_affected || child->affected!=2) {
	  if (miss_f) {
	    if (fill_in(child->markers+ia1, f->father+ia1, f->mother+ia1))
	      return 1;
	  }
	  else {
	    if (fill_in(child->markers+ia1, f->mother+ia1, f->father+ia1))
	      return 1;
	  }
	}
      }
    }
  }
  return 0;
}

/* 
   Expand family into multiple copies, one copy for each offspring affected 
   These will be consecutive as the extra records are inserted in the list
   between f and f->next. Returns pointer to the last family created.
*/

Family *expand_family(Family *f, int m) {
  int naff, iaff, jaff;
  Offspring *child;
  Family *f_old, *f_new, *f_end;
  /* First find out how many affected offspring there are */
  for (naff=0, child = f->children; child; child=child->next) {
    if (child->affected==2) naff++;
  }
  /* Now make (naff -1) copies of the family */
  f_end = f->next;
  f_old = f;
  for(iaff=1; iaff<naff; iaff++) {
    f_new = copy_family(f, m);
    if (!f_new) {
      warn("Not enough memory to copy family", f);
      f->next = f_end;
      return 0;
    }
    f_old->next = f_new;
    f_old = f_new;
  }
  /* Now, in every duplicate family, reset all affected bar 1 to 0 */
  if (naff>1) {
    for (f_new=f, iaff=0; f_new; f_new=f_new->next, iaff++) {
      for (jaff=0, child=f_new->children; child; child=child->next) {
	if (child->affected==2) {
	  if (jaff!=iaff)	child->affected = 0;
	  jaff++;
	}
      }
    }
  }
  /* Rejoin last element to rest of list */
  f_old->next = f_end;
  return f_old;
}

/* Create a family record */

Family *new_family(int m) {
  int locus, twom, ia;
  Family *res;
  /* Storage allocation */
  res = (Family *) malloc(sizeof(Family));
  if (!res) return (Family *) 0;
  res->father = (int *) calloc(2*m, sizeof(int));
  if (!res->father) return (Family *) 0;
  res->mother = (int *) calloc(2*m, sizeof(int));
  if (!res->mother) return (Family *) 0;
  res->phase_f = (int *) calloc(m, sizeof(int));
  if (!res->phase_f) return (Family *) 0;
  res->phase_m = (int *) calloc(m, sizeof(int));
  if (!res->phase_m) return (Family *) 0;
  /* Initialization */
  twom = 2*m;
  res->check = 0;
  res->pedigree = res->father_id = res->mother_id = 0;
  for (locus = 0; locus<m; locus++) 
    res->phase_f[locus] = res->phase_m[locus] = 0;
  for (ia=0; ia<twom; ia++) 
    res->father[ia] = res->mother[ia] = 0;
  res->children = (Offspring *) 0;
  res->next = (Family *) 0;
  return res;
}

/* Set affected status for all affected offspring other than first to zero */

void use_only_first(Family *f) {
  Offspring *child;
  int reset;
  for (reset=0, child=f->children; child; child=child->next) {
    if (child->affected==2) { 
      if (reset) { 
	child->affected = 0;
      }
      else {
	reset = 1;
      }
    }
  }
}

/* Copy family record --- except for phase vectors */

Family *copy_family(Family *f, int m) {
  int i, twom;
  Family *res;
  Offspring *child, *twin, *previous;
  twom = 2*m;
  res = new_family(m);
  if (!res) return (Family *) 0;
  res->check = f->check;
  res->pedigree = f->pedigree;
  res->father_id = f->father_id;
  res->mother_id = f->mother_id;
  for (i=0; i<twom; i++) {
    res->father[i] = f->father[i];
    res->mother[i] = f->mother[i];
  }
  /* Now copy children */
  for (previous=(Offspring *) 0, child=f->children; child; child=child->next){
    twin = copy_child(child, m);
    if (!twin) return (Family *) 0;
    if (previous) 
      previous->next = twin;
    else 
      res->children = twin;
    previous = twin;
  }
  return res;
}

/* Delete family record */

void del_family(Family *f) {
  Offspring *child, *chdel;
  free(f->father);
  free(f->mother);
  free(f->phase_f);
  free(f->phase_m);
  child = f->children;
  while (child) {
    chdel = child;
    child = child->next;
    del_child(chdel);
  }
}

/* Create new child record */

Offspring *new_child(int m) {
  int locus, twom, ia;
  Offspring *res;
  /* Storage allocation */
  res = (Offspring *) malloc(sizeof(Offspring));
  if (!res) return (Offspring *) 0;
  res->markers = (int *) calloc(2*m, sizeof(int));
  if (!res->markers) return (Offspring *) 0;
  res->ivec_f = (int *) calloc(m, sizeof(int));
  if (!res->ivec_f) return (Offspring *) 0;
  res->ivec_m = (int *) calloc(m, sizeof(int));
  if (!res->ivec_m) return (Offspring *) 0;
  /* Initialization */
  twom = 2*m;
  for (locus=0; locus<m; locus++) 
    res->ivec_f[locus] = res->ivec_m[locus] = 0;
  for (ia=0; ia<twom; ia++) 
    res->markers[ia] = 0;
  res->id = res->affected = res->f_tr = res->m_tr = res->sex = 0;
  res->qt = 0.0;
  res->next = (Offspring *) 0;
  return res;
}

/* Copy a child record --- except for transmission and inheritance vectors */

Offspring *copy_child(Offspring *child, int m) {
  int ia, twom;
  Offspring *res;
  res = new_child(m);
  if (!res) return (Offspring *) 0;
  res->id = child->id;
  res->affected = child->affected;
  res->sex = child->sex;
  res->qt = child->qt;
  twom = 2*m;
  for (ia=0; ia<twom; ia++) {
    res->markers[ia] = child->markers[ia];
  }
  return res;
}

/* Delete a child record */

void del_child(Offspring *child) {
  free(child->markers);
  free(child->ivec_f);
  free(child->ivec_m);
  free(child);
}

/* Print a one-line summary of the family to a file */

void show_family(Family *f) {
  Offspring *child;
  if (f) {
    REprintf(" %d: %d + %d / ", 
	    f->pedigree, f->father_id, f->mother_id);
    for (child=f->children; child; child=child->next) {
      REprintf(" %d", child->id);
      if (child->affected==2) 
	REprintf("*");
      if (child->next) 
	REprintf(",");
    }
    REprintf("\n");
  }
  else {
    REprintf("*** empty family ***\n");
  }
}

/* Print the family data in full */

void print_family(Family *f, int m, FILE *stream) {
  int i, j, ia1, ia2;
  Offspring *child;
  if (!f) return;
  fprintf(stream, "Pedigree %8d:\n     Father      Mother", f->pedigree);
  for (i=0, child=f->children; child && i<4; child=child->next, i++)
    fprintf(stream, "  Offspring%3d", i+1);
  if (child)
    fprintf(stream, ">\n");
  else
    fprintf(stream, "\n"); 
  fprintf(stream, "%8d Ph %8d Ph", f->father_id, f->mother_id);
  for (i=0, child=f->children; child && i<4; child=child->next, i++) {
    if (child->affected==2) 
      fprintf(stream, "   %8d*Iv", child->id);
    else
      fprintf(stream, "   %8d Iv", child->id);
  }
  fprintf(stream, "\n");
  /* Print markers, plus phase and inheritance vectors */
  for (ia2=1, ia1=j=0; j<m; j++, ia1+=2, ia2+=2) {
     fprintf(stream, "%5d%4d%2d%6d%4d%2d",
	    f->father[ia1], f->father[ia2], f->phase_f[j],
	    f->mother[ia1], f->mother[ia2], f->phase_m[j]);
    for (i=0, child=f->children; child && i<4; child=child->next, i++) {
      fprintf(stream, "%7d%4d%2d%1d", 
	      child->markers[ia1], child->markers[ia2],
	      child->ivec_f[j], child->ivec_m[j] );
    }
    fprintf(stream, "\n");
  }
  /* Print haplotype transmission */
  fprintf(stream,"Haplotype transmission:");
  for (i=0, child=f->children; child && i<4; child=child->next, i++) 
    fprintf(stream, "            %d%d", child->f_tr, child->m_tr);
  fprintf(stream, "\n");
}

/* Print a warning to stderr */

void warn(char *message, Family *f) {
  REprintf(message);
  REprintf(": ");
  show_family(f);
}

/* Count families */

int count_families(Family *first) {
  int res;
  Family *f;
  for (res=0, f=first; f; f=f->next, res++){}
  return res;
}

/* Count offspring */

int count_offspring(Family *first, int affected_only) {
  int res;
  Family *f;
  Offspring *child;
  for (res=0, f=first; f; f=f->next) {
    for (child=f->children; child; child=child->next) {
      if (affected_only) {
	if (child->affected==2) res++;
      }
      else {
	res++;
      }
    }
  }
  return res;
}

/* ============================ Local functions =============================*/


/* 
   Which allele is inherited from parent. Returns:
   0:  Neither
   1:  First
   2:  Second
   3:  Either
*/

int inherited(int allele, int p_gtype[2]) {
  int res;
  if (allele==0) return 3;
  res = 0;
  if (*p_gtype == 0 || *p_gtype == allele) res++;
  if (*(p_gtype+1) == 0 || *(p_gtype+1) == allele) res+=2;
  return res;
}

/*
  Is transmission pattern consistent with offspring genotype?
*/

int poss_tr(int allele_f, int allele_m, int o_gtype[2]) {
  int o1, o2, f1, f2, m1, m2, res;
  o1 = o_gtype[0];
  o2 = o_gtype[1];
  f1 = !allele_f || !o1 || (allele_f==o1);
  f2 = !allele_f || !o2 || (allele_f==o2);
  m1 = !allele_m || !o1 || (allele_m==o1);
  m2 = !allele_m || !o2 || (allele_m==o2);
  res = (f1+f2) && (m1+m2) && (f1+m1) && (f2+m2);
  return res;
}

/*
  Given a parental phase vector and an inheritance vector, determines the 
  transmission. Returns
  0: if transmission is ambiguous
  1: if in-phase haplotype transmitted,
  2: if out-of-phase haplotype transmitted
  -L:if recombination at locus L
  If 1 or 2 is returned and parental phase vector contains ambiguities 
  (zero entries), then complete them where possible.
  If phase vector is empty and inheritance vector is not, returns 1 and sets
  phase vector = inheritance vector.
*/

int trans(int phase[], int inhvec[], int m) {
  int locus, in_ph=0, resolved, nph, nin;
  /* Find loci that can resolve phase */
  for (nph=nin=resolved=locus=0; locus<m; locus++) {
    if (phase[locus]) nph++;
    if (inhvec[locus]) nin++;
    if (phase[locus] && inhvec[locus]) {
      /* If already resolved, check for recombination */
      if (resolved) {
	if (in_ph && (phase[locus]!=inhvec[locus])) return -(1+locus); 
      }
      /* else resolve phase */
      else {
	in_ph = (phase[locus]==inhvec[locus]);
	resolved = 1;
      }
    }
  }
  /* If empty phase vector and non-empty inheritance vector, 
     resolve as in-phase */
  if (nph==0 && nin>0) {
    for (locus=0; locus<m; locus++) {
      phase[locus] = inhvec[locus];
    }
    return 1;
  }
  if (!resolved) return 0;
  /* Fill in any gaps in the phase vector */
  for (locus=0; locus<m; locus++) {
    if (inhvec[locus] && !phase[locus]) {
      phase[locus] = in_ph? inhvec[locus]: 3 - inhvec[locus];
    }
  }
  if (in_ph) return 1;
  return 2;
}

/* 
   Fill in alleles in unknown parental genotype given child genotype and
   known parental genotype. Returns 1 if error/ex-paternity
*/

int fill_in(int child[2], int unknown[2], int known[2]) {
  int which;
  /* Is first allele present and not inherited from known genotype? */
  if (*child && !inherited(*child, known)) {
    which = 0;
    /* Go on to put child[0] into first available slot */
  }
  else {
    which = 1;
    /* Now child[1] is the candidate for imputation */
  }
  /* Is second allele present and not inherited from known genotype? */
  if (*(child+1) && !inherited(*(child+1), known)) {
    if (which==0) { /* Both alleles can't come from other parent! */
       return 1;
    }
    /* Go on to put put child[1] into first available slot */
  }
  else {
    /* If neither qualifies and heterozygous, then no action */
    if (which==1 && (*child != *(child+1))) return 0;
  }
  /* Put child[which] in first available slot */
  which = *(child+which);
  if (!*unknown) {
    *unknown = which;
  }
  else if (*unknown != which) {
    if (!*(unknown+1)) {
      *(unknown+1) = which;
    }
    else if (*(unknown+1) != which) { 
      return 1;
    }
  }
  return 0;
}
