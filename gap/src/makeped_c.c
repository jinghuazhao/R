/*****************************************************************************/
/*                                                                           */
/*  Program: MAKEPED.C                                                       */
/*                                                                           */
/*                                                                           */
/*  Used to convert pedigree files missing sib pointers to pedigree files    */
/*  with sib pointers.                                                       */
/*                                                                           */
/*  Revision history:                                                        */
/*                                                                           */
/*  6-23-88  The ind_lookup routine stopped looking for an id too soon.      */
/*           This showed up when strings where used as id's and the pedigree */
/*           was entered in the reverse order from the way you would draw    */
/*           a pedigree, that is, an individual drawn at bottom of a pedigree*/
/*           was listed as the first individual in the data file.            */
/*                                                                           */
/* 10-05-88  Character spacing of output when loops created an new individual*/
/*           has been corrected.                                             */
/*                                                                           */
/* 10-09-88  Loops have been corrected so that offspring refer to the newly  */
/*           created individual as a parent instead of the original parent.  */
/*                                                                           */
/* 24-12-88  Loops have been corrected so that nextoffspring points are      */
/*           kept for the original individual instead of new individual.     */
/*                                                                           */
/*  8-30-89  Removed a lot of code dealing with ped ids and person ids as    */
/*           integers.                                                       */
/*           Ids are read in as strings; original pedigree ids are only      */
/*           used if ALL pedigree ids consist of digits.                     */
/*           Allows for multiple loops, with checks for: if a loop person is */
/*           the proband, then the proband must be in the first loop.        */
/*           Loop_file and Proband_file names are read with gets() instead   */
/*           of fscanf().                                                    */
/*  9-14-90  Changed Version number to 2.21 and re-compiled with large       */
/*           memory model.  Deleted onscreen references to my name (Daniel   */
/*           E. Weeks).  Change pifile and pofile to s_byte from u_byte.     */
/*  1-16-92  Added a feature to check duplicated individual IDs which are    */
/*           supposed to be unique.  Added by Xiaoli Xie.                    */
/*  1-19-04  No further warnings for R port. Jing hua Zhao                   */
/*****************************************************************************/
#define MAKEPED_VERSION 2.21
#include <R.h>
#include <stdio.h>
#include <ctype.h>
/* #include <Storage.h> */
#include <stdlib.h>  /* Macintosh addition */
/* #include <unix.h> */  /* Microsoft C doesn't like this */
#include <string.h>
#ifdef MICROSOFT_C
/*  #include <malloc.h> malloc/malloc.h for Mac OS X 10.3 (Panther)*/
#endif

#ifdef TURBO_C
  #include <alloc.h>
#endif

#define	FALSE			0
#define	TRUE			1

typedef unsigned char   u_byte;
typedef          char   s_byte;
typedef unsigned short  u_word;
typedef          short  s_word;
typedef          long   s_long;
typedef unsigned int    u_intg;
typedef          int    s_intg;
typedef          float  s_real;
typedef          double d_real;

        /* User defined constants */

#define maxallchars  200      /* maximum chars used in phenotypic data       */
#define maxind       8001     /* maximum number of individuals               */
#define maxped       2401      /* maximum number of pedigrees                 */
#define maxname      11       /* maximum number of chars in ped or person id */
#define maxloop      4        /* The total number of loops over all pedigr.  */
                              /* cannot exceed maxped*maxloop                */
char BELL='\007';

        /* File pointers */

FILE *pedfile;                /* pointer to pedigree information             */
FILE *pedout;                 /* pointer to new pedigree output file         */

        /* Others */

#define max_filespec   80     /* max number of chars in a filename           */

struct phenotype {
  u_byte phen_chars[maxallchars];
};

struct ind{
  s_byte oldped_s[maxname]; /* original pedigree read from pedfile if string */      
  s_intg oldped;          /* original pedigree read from pedfile if integer  */      
  s_byte oldid_s[maxname];/* original person id read from pedfile if string  */ 
  s_intg oldid;           /* original person id read from pedfile if integer */
  s_intg ped;             /* new pedigree id                                 */
  s_intg id;              /* new person id                                   */
  s_intg paid;            /* new paternal id                                 */
  s_intg maid;            /* new maternal id                                 */
  s_intg offid;           /* new first offspring id                          */
  s_intg npaid;           /* new next paternal offspring id                  */
  s_intg nmaid;           /* new next maternal offspring id                  */
  s_intg sex;             /* sex of person                                   */
  s_intg proband;         /* proband id = 1; loop id >= 2; can be -1 (see below) */
  struct ind *pa;         /* pointer to fathers ind structure                */
  struct ind *ma;         /* pointer to mothers ind structure                */
  struct ind *foff;       /* pointer to first offsprings ind structure       */
  struct ind *nextpa;     /* pointer to next paternal offspring ind struct   */
  struct ind *nextma;     /* pointer to next maternal offspring ind struct   */
  s_intg generations;     /* number of generations below this individual     */
  struct phenotype *phen; /* pointer to persons phenotype structure          */
  s_intg male;            /* boolean flag for sex of individual              */
  u_intg is_parent;       /* flag used to find people with no family         */
} ;

        /* arrays of individuals */
struct ind *person[ maxind ];  /* array of pointers to ind structures        */
struct ind *proband[ maxped ]; /* array of pointers to pedigree probands     */

        /* others */

s_intg nuped;
s_intg next_id;             /* Person id's used if originals were strings.  */
s_intg nuperson;            /* Number of people read so far in a pedigree.  */
s_intg totperson;           /* Total num persons in pedigree file.          */
s_intg probands[maxped];    /* Indexes of people assigned probands.         */
s_intg loops[maxped*maxloop];       /* Indexes of people assigned loops.            */
u_intg found_error;         /* Error flag for sex and phenotype check.      */
s_byte pifile[max_filespec];/* Pedigree input file.                         */
s_byte pofile[max_filespec];/* Revised pedigree output file.                */
s_intg biggest_p_id;        /* Largest pedigree id found or created.        */
s_intg biggest_i_id;        /* Largest individual id found or created.      */
s_intg loop_i = 2;          /* Counter for loop ids in proband field */

s_byte curped_s[maxname];   /* Current pedigree id                          */
s_intg lineperson;	    /* Number of people per family read so far      */ 
s_byte lineind[maxind][maxname]; /*Array stores the people of current family*/

 /* The following are set to true if the pedigree has already been cleared  */
 /* of probands set by auto_probands.  You only want to clear probands from */
 /* a pedigree once so that a user can set more than one proband/pedigree   */
 /* if he wants to.                                                         */

u_byte cleared[maxped] = {FALSE};

 /* ped_integers is true only if all pedigree ids consist of digits; only */
 /* in this case, will the original pedigree ids used in the output file. */

u_byte ped_integers;

/****************************************************************************/
/*                                                                          */
/*			   read_pedigree   				    */
/*                                                                          */
/* Prompts for the identifier of a pedigree                                 */
/*                                                                          */
/****************************************************************************/

void read_pedigree(pedigree_s)
 s_byte pedigree_s[];   /* string id's  */
{
  Rprintf("\n\tPedigree   -> ");
  fscanf(stdin,"%s",pedigree_s);
}
/****************************************************************************/
/*                                                                          */
/*			   read_person   				    */
/*                                                                          */
/* Prompts for the identifier of a person                                   */
/*                                                                          */
/****************************************************************************/

void read_person(person_s)
 s_byte person_s[];   /* string id's  */
{
  Rprintf("\tPerson     -> ");
  fscanf(stdin,"%s",person_s);
}
/****************************************************************************/
/*                                                                          */
/*                            ind_lookup                                    */
/*                                                                          */
/****************************************************************************/

s_intg ind_lookup(name,sequence)
  s_byte name[];
  s_intg sequence;  /* the number of people in previous pedigrees */
{

/* Search through all people in this pedigree that have been read in so far */
/* and find the one whose original id matches that in the name parameter.   */
/* If one is found then return his id else return 0.                        */

s_intg i;

i = 1;

while (i <= nuperson) {
  if (person[sequence+i] == NULL)  return(0);
  else
   if (!strcmp(person[sequence + i]->oldid_s,name) ) /* if they match */
        return(person[sequence + i]->id);
  else
    i++;
}
return(0);
}

/****************************************************************************/
/*                                                                          */
/*                            chk_dupli                                     */
/*                                                                          */
/****************************************************************************/

s_intg chk_dupli(name)
  s_byte name[];
/*  s_intg sequence;   the number of people in previous pedigrees */
{

/* Search through all people in this pedigree that have been read in so far */
/* and see if this individual id is unique. If it is unique return 0 else   */
/* return 1.                                                                */

s_intg i;

i = 1;
while (i <= lineperson) {
  if (lineind[i][0] == '\0')  return(0);
  else
   if (!strcmp(lineind[i],name) ) /* if they match */
    { 
      Rprintf("\nWARNING! Individual id. %s in family %s is duplicated%c\n",name,curped_s,BELL);
      return(1);
     }
  else
    i++;
}
return(0);
}


/*****************************************************************************/
/*                                                                           */
/*                             getind                                        */
/*                                                                           */
/*****************************************************************************/

void getind(id,sequence,newped_s,nuped)
     s_intg *id;
     s_intg sequence;
     s_byte newped_s[];
     s_intg nuped;

{
/*s_intg temp_i;                  Holds id if id's are integers. */
  s_byte temp_s[maxname];      /* Holds id if id's are strings.  */
  s_intg found_id;
  s_intg dupli_id;             /* added by Xiaoli Xie            */

  /* Read the persons id and convert it to an initial index. */

    fscanf(pedfile,"%s",temp_s);
    dupli_id = chk_dupli(temp_s);         /* yes:1 no:0 */
    if (!dupli_id) {
    strcpy(lineind[lineperson],temp_s);
    lineperson++;
    }
    if (temp_s[0]=='0' && temp_s[1]=='\0') *id = 0;
    else {
      found_id = ind_lookup(temp_s,sequence);
      if (found_id) *id = found_id;
      else *id = next_id;
    }
  
/* Should have checking if !ind_integers
    if ((ind_integers) && (temp_i > maxind)) {
     REprintf("\nERROR: Ped: %d  Per: %d - ",newped_i,temp_i);
     REprintf("maximum id of %d exceeded.\n",maxind - 1);
     error("%d",1);
    }  */


  /* If the initial index is not zero then compute the final index, */
  /* and allocate memory for their record if not already done, and  */
  /* initialize their record fields.                                */

  if  (*id != 0) {
      *id += sequence;
      if (person[*id] == NULL) {
	person[*id] = (struct ind *) calloc ( 1, sizeof( struct ind ) );

	if ( person[*id] == NULL ) {
	  error("\nERROR: Cannot allocate memory.\n");
	}

        strcpy(person[*id]->oldped_s,newped_s);
        strcpy(person[*id]->oldid_s,temp_s);
		person[*id]->id       = next_id;
        if (next_id > biggest_i_id) biggest_i_id = next_id;
        next_id++;
        if (nuped > biggest_p_id) biggest_p_id = nuped;

      person[*id]->ped        = nuped;
      person[*id]->paid       = 0;
      person[*id]->maid       = 0;
      person[*id]->offid      = 0;
      person[*id]->npaid      = 0;
      person[*id]->nmaid      = 0;
      person[*id]->pa         = NULL;
      person[*id]->ma         = NULL;
      person[*id]->foff       = NULL;
      person[*id]->nextpa     = NULL;
      person[*id]->nextma     = NULL;
      person[*id]->proband    = 0; /* Initialize proband field */
      nuperson++;

      }
    }
}


/*****************************************************************************/
/*                                                                           */
/*                             getindpa                                      */
/*                                                                           */
/*****************************************************************************/

void getindpa(id,sequence,newped_s,nuped)
     s_intg *id;
     s_intg sequence;
     s_byte newped_s[];
     s_intg nuped;

{
/*s_intg temp_i;                  Holds id if id's are integers. */
  s_byte temp_s[maxname];      /* Holds id if id's are strings.  */
  s_intg found_id;

  /* Read the persons id and convert it to an initial index. */

    fscanf(pedfile,"%s",temp_s);
    if (temp_s[0]=='0' && temp_s[1]=='\0') *id = 0;
    else {
      found_id = ind_lookup(temp_s,sequence);
      if (found_id) *id = found_id;
      else *id = next_id;
    }
  
/* Should have checking if !ind_integers
    if ((ind_integers) && (temp_i > maxind)) {
     REprintf("\nERROR: Ped: %d  Per: %d - ",newped_i,temp_i);
     REprintf("maximum id of %d exceeded.\n",maxind - 1);
     error("%d",1);
    }  */


  /* If the initial index is not zero then compute the final index, */
  /* and allocate memory for their record if not already done, and  */
  /* initialize their record fields.                                */

  if  (*id != 0) {
      *id += sequence;
      if (person[*id] == NULL) {
	person[*id] = (struct ind *) calloc ( 1, sizeof( struct ind ) );

	if ( person[*id] == NULL ) {
	  error("\nERROR: Cannot allocate memory.\n");
	}

        strcpy(person[*id]->oldped_s,newped_s);
        strcpy(person[*id]->oldid_s,temp_s);
		person[*id]->id       = next_id;
        if (next_id > biggest_i_id) biggest_i_id = next_id;
        next_id++;
        if (nuped > biggest_p_id) biggest_p_id = nuped;

      person[*id]->ped        = nuped;
      person[*id]->paid       = 0;
      person[*id]->maid       = 0;
      person[*id]->offid      = 0;
      person[*id]->npaid      = 0;
      person[*id]->nmaid      = 0;
      person[*id]->pa         = NULL;
      person[*id]->ma         = NULL;
      person[*id]->foff       = NULL;
      person[*id]->nextpa     = NULL;
      person[*id]->nextma     = NULL;
      person[*id]->proband    = 0; /* Initialize proband field */
      nuperson++;

      }
    }
}

/*****************************************************************************/
/*                                                                           */
/*                             getphenotype                                  */
/*                                                                           */
/* This version of getphenotype simply reads in the phenotypic data as a     */
/* string and assigns it to the persons phen->phen_chars field.  This should */
/* be used when no interpretation of phenotypic data is intended.            */
/*                                                                           */
/*****************************************************************************/

void getphenotype(id)
     s_intg *id;

{
  s_intg i;
  s_byte c;

  person[*id]->phen =
    (struct phenotype *) calloc ( 1, sizeof( struct phenotype ) );
  if ( person[*id]->phen == NULL ) {
      error("\nERROR: Cannot allocate memory.\n");
  }

  i = 0;
  c = getc(pedfile);
  while(( c != '\n') && (!feof(pedfile))) {
    person[*id]->phen->phen_chars[i++] = c;
    c = getc(pedfile);
  }

  person[*id]->phen->phen_chars[i] = '\0';

}
/*****************************************************************************/
/*                                                                           */
/*                             readped                                       */
/*                                                                           */
/*****************************************************************************/

void readped()
{
  s_intg i;
/*s_intg j;*/
  s_intg newid;
  s_intg sex;
/*s_intg profield;*/
  s_intg sequence;                /* number of people read in pedigree */
  s_intg thisone;
  s_byte thisped_s[maxname];      /* string pedigree id  */
  s_byte newped_s[maxname];       /* string pedigree id  */
/*s_byte c;*/

  /* initialize some stuff */

  totperson      = 0;
  sequence       = 0;
  nuperson       = 0;
  nuped          = 1;
  next_id        = 1;
  lineperson =1; 	 /*added by xie*/
  proband[nuped] = NULL;
  for( i=0; i<maxind; i++) person[i] = NULL;
  for( i=0; i<maxind; i++) lineind[i][0] ='\0';  /*   */

    /* Read first pedigree number*/

  rewind(pedfile);
  fscanf(pedfile,"%s",newped_s);
  strcpy(thisped_s,newped_s);
  strcpy(curped_s,newped_s);
  while (! feof(pedfile)) {

    /* Get person. */

    getind(&thisone,sequence,thisped_s,nuped); 

    /* Get persons father. */

    getindpa(&newid,sequence,thisped_s,nuped);
    /* Note: the single = sign does ASSIGNMENT of person[newid] to pa */
    if ((person[thisone]->pa = person[newid]) != NULL)
      person[thisone]->paid   = person[newid]->id;

    /* Get persons mother. */

    getindpa(&newid,sequence,thisped_s,nuped);
    if ((person[thisone]->ma = person[newid]) != NULL)
      person[thisone]->maid   = person[newid]->id;

    fscanf(pedfile,"%d",&sex);
    person[thisone]->sex = sex;
    if (sex == 1) person[thisone]->male = TRUE;
    else person[thisone]->male = FALSE;

    getphenotype(&thisone);

    /* Read in the next the pedigree id, if it's the start of the next */
    /* pedigree then update a few things and keep going.               */

    if (!feof(pedfile)) {
        fscanf(pedfile,"%s",newped_s);
        if (strcmp(thisped_s,newped_s) != 0) {
        sequence += nuperson;
        nuperson = 0;
        lineperson = 1; /* added by xie*/
        for(i=0;i<maxind;i++) lineind[i][0]='\0';
        strcpy(thisped_s,newped_s);
        strcpy(curped_s,newped_s);
        nuped++;
        next_id = 1;  /* next id used if person id's were strings */
        }
      
    }
  }
  totperson = nuperson + sequence;
  /* Note: C uses zero-based arrays, so the indices go from 0 to maxind-1, */
  /* but we are not using person zero. */
  if (totperson > maxind-1) {
   error("\nERROR: maximum number %d of individuals exceeded \n",maxind-1);
  }
}

/* *******************************************************************  */
void pointers()
{
  struct ind *q;
  s_intg      i;
  s_intg      count;      /* number of people in current pedigree       */
  s_intg      ped_count;  /* number of people in all previous pedigrees */
  s_intg      ped_id;     /* current pedigree number                    */

     /* Note: these variable are easy to confuse... */
     /*                                             */
     /*           id's    pointer's                 */
     /*           ------  ---------                 */
     /*           paid    pa                        */
     /*           maid    ma                        */
     /*           offid   foff                      */
     /*           npaid   nextpa                    */
     /*           nmaid   nextma                    */

  count      = 0;
  ped_count  = 0;
  ped_id     = 0;

  for(i=1; i<=totperson; i++){
   if (person[i] != NULL) {

     if (person[i]->ped != ped_id) {  /* if beginning of next pedigree... */
      ped_id    = person[i]->ped;
      ped_count += count;
      count     = 0;
    }
     ++count;

     if ( person[i]->paid != 0 ) {
       q = person[ person[i]->paid + ped_count ]; /* get pointer to father */

     /* If the fathers first offspring is not set then set it to person i */

       if (q->offid == 0) {
	 q->offid = i - ped_count;
	 q->foff  = person[i];
       }
       else {

     /* If the fathers first offspring is set then work your way down the */
     /* next paternal chain, starting from the first offspring, until it  */
     /* runs out. Then set the next paternal of the last person in the    */
     /* chain to the original person i that you started out with.         */

	 q = person[q->offid + ped_count];
	 while (q->nextpa != NULL) q = person[q->npaid + ped_count];
         q->npaid  =  i - ped_count;
	 q->nextpa = person[i];
       }
     }

     /* Repeat the above procedure for the mothers side. */

     if ( person[i]->maid != 0 ) {
       q = person[ person[i]->maid + ped_count ]; /* get pointer to mother */
       if (q->offid == 0) {
	 q->offid =  i - ped_count;
	 q->foff  = person[i];
       }
       else {
	 q = person[q->offid + ped_count];
	 while (q->nextma != NULL) q = person[q->nmaid + ped_count];
         q->nmaid  =  i - ped_count;
	 q->nextma = person[i];

       }
     }
   }
 }
}
/*****************************************************************************/
/*                                                                           */
/*                          save_loops                                       */
/*                                                                           */
/*****************************************************************************/

void save_loops(count)
    s_intg count;
{
  s_intg i;
  s_byte response;
  s_byte loop_file[max_filespec];
  FILE   *loopf;


  Rprintf("\n\nDo you want these selections saved ");
  Rprintf("for later use?  (y/n) -> ");
  fscanf(stdin,"%1s",&response);

  if ((response == 'y') || (response == 'Y')) {
    loop_file[0] = '\0';
    Rprintf("\nEnter filename -> ");
    while ( loop_file[0] == '\0' ) {
      fgets(loop_file,max_filespec,stdin);
    }
    if ( (loopf = fopen(loop_file,"w")) == NULL) {
      REprintf("\nERROR: Cannot open file %s\n",loop_file);
    }
    else {
      for(i=0; i<count; i++) {
         fprintf(loopf,"%s ",person[ loops[i] ]->oldped_s); 
         fprintf(loopf,"%s\n",person[ loops[i] ]->oldid_s);
      }
      fclose(loopf);
    }
  }
 }
/****************************************************************************/
/*                                                                          */
/*                            largest_id                                    */
/*                                                                          */
/* Given the index of a person, this routine returns the largest integer id */
/* found in that persons pedigree.                                          */
/*                                                                          */
/****************************************************************************/

s_intg largest_id(person_index)
   s_intg person_index;
{

/* s_intg pedigree_number;*/
   s_intg largest;
   s_intg i;

   largest = person[person_index]->id;

   i = person_index -1;
   while (( i >= 1 ) && (person[i]->ped == person[person_index]->ped)) {
     if (person[i]->id > largest )
       largest = person[i]->id;
     --i;
   }

   i = person_index +1;
   while((person[i] != NULL) && (person[i]->ped == person[person_index]->ped)){
     if (person[i]->id > largest )
       largest = person[i]->id;
     ++i;
   }

   return(largest);
 }
/*****************************************************************************/
/*                                                                           */
/*                            add_loop                                       */
/*                                                                           */
/*****************************************************************************/

void add_loop(start_of_ped,old)
  s_intg start_of_ped;
  s_intg old;
{
  s_intg i,max_i;
  s_intg new;
  s_intg next_possible_id;
  s_intg pedigree;

  /* Scan through the current pedigree to find the maximum loop_i  */
  /* assigned so far.  Increment this by one.  If none have been   */
  /* assigned, start with loop_i = 2. Do this before the new space */ 
  /* has been opened up                                            */

  max_i = 1;
  i = start_of_ped;
  pedigree = person[start_of_ped]->ped;
  while ((i <= totperson) && (pedigree == person[i]->ped)) {
   if(person[i]->proband > max_i)
     max_i = person[i]->proband;
   i++;
  }
  loop_i = max_i + 1;
  
  /* Get next possible id for this pedigree */

  next_possible_id = largest_id(old) + 1;
  if (next_possible_id > biggest_i_id) biggest_i_id = next_possible_id;

  /* Open a slot in the person array to insert a new person. */

  i = totperson;
  while(i > old) {
   person[i+1] = person[i];
   i--;
  }
  new = i + 1;
  totperson++;
  if (totperson > maxind-1) {
   error("\nERROR: maximum number %d of individuals exceeded \n",maxind-1);
  }

  person[new] = (struct ind *) calloc ( 1, sizeof( struct ind ) );
  if ( person[new] == NULL ) {
    error("\nERROR: Cannot allocate memory.\n");
  }

  /* Copy the original record to the new record  */
  /* but 0 out the parents of the new record and */
  /* and the children of the old record.         */

    strcpy(person[new]->oldped_s,person[old]->oldped_s);
    strcpy(person[new]->oldid_s,person[old]->oldid_s);

  person[new]->id         = next_possible_id;
  person[new]->ped        = person[old]->ped;
  person[new]->paid       = 0;
  person[new]->maid       = 0;
  person[new]->pa         = NULL;
  person[new]->ma         = NULL;
  person[new]->offid      = person[old]->offid;
  person[new]->foff       = person[old]->foff;
  person[new]->nextpa     = NULL;
  person[new]->nextma     = NULL;
  person[new]->sex        = person[old]->sex;
  person[new]->proband    = loop_i; 
  person[new]->phen       = person[old]->phen;

  /* zero out the children of the original record */
  /* and assign the proband.                      */

  person[old]->offid      = 0;
  person[old]->npaid      = 0;
  person[old]->nmaid      = 0;
  person[old]->foff       = NULL;
  person[old]->proband    = loop_i; 

  /* Scan through entire pedigree and look for people who had a parent */
  /* that was the original and switch them to the newly created parent */

       i = start_of_ped;
       pedigree = person[start_of_ped]->ped;

       while( (i<=totperson) &&
            (pedigree == person[i]->ped)){

           if (person[i]->paid == person[old]->id) {
              person[i]->paid = person[new]->id;
              person[i]->pa = person[new];
            }
	   if (person[i]->maid == person[old]->id) {
              person[i]->maid = person[new]->id;
              person[i]->ma = person[new];
           }
         i++;
        } /* end of while */
}

/*****************************************************************************/
/*                                                                           */
/*                             file_loops                                    */
/*                                                                           */
/*****************************************************************************/
#ifdef executable
void file_loops()
#else
void file_loops(char **loopfile)
#endif
{
/*s_byte response;*/
/*s_intg pedigree;*/
/*s_intg pedigree_i;*/
  s_byte pedigree_s[maxname];
/*s_intg person_i;*/
  s_byte person_s[maxname];
  s_intg i;
/*s_intg from_file;*/
  s_intg found_start_of_ped;
  s_intg start_of_ped=0;
  s_intg found_person;
/*s_intg good_format;*/
/*s_byte c;*/
  FILE   *loopf;
#ifdef executable 
  s_byte loop_file[max_filespec];
  loop_file[0] = '\0';
  Rprintf("\nEnter filename -> ");
  while ( loop_file[0] == '\0' ) {
    fgets(loop_file,max_filespec,stdin);
  }
  if ( (loopf = fopen(loop_file,"r")) == NULL) {
   error("\nERROR: Cannot open file %s\n",loop_file);
  }
#else
   if ((loopf = fopen(*loopfile,"r")) == NULL) {
       error("\nERROR: Cannot open file %s\n",*loopfile);
   }
#endif

  while(!feof(loopf)) {

   fscanf(loopf,"%s",pedigree_s);
   fscanf(loopf,"%s",person_s);

   if (!feof(loopf)) {

  /* convert to integer if needed */

  found_person = FALSE;
  found_start_of_ped = FALSE;
  i = 1;
  while(( i <= totperson) && (!found_person)) {

    if ((!found_start_of_ped) && (strcmp(pedigree_s,person[i]->oldped_s)==0)) {
      start_of_ped = i;
      found_start_of_ped = TRUE;
    }
    if ( (strcmp(pedigree_s,person[i]->oldped_s) == 0) &&
         (strcmp(person_s,person[i]->oldid_s ) == 0)) {
           found_person = TRUE;
           add_loop(start_of_ped,i);
    }
    i++;
	if (( i>totperson) && (!found_person)) {
           error("\nERROR: Ped: %s Per: %s - Not found, check loop file.\n",
		  pedigree_s,
		  person_s);
	}
      }
     }
    }
     fclose(loopf);
}

/*****************************************************************************/
/*                                                                           */
/*                             some_loops                                    */
/*                                                                           */
/*****************************************************************************/

void some_loops()
{
/*s_byte response;*/
/*s_intg person_i;*/
  s_byte person_s[maxname];
/*s_intg pedigree_i;*/
  s_byte pedigree_s[maxname];
  s_intg pedigree=0;
  s_intg i,/*j,k,*/ii;
  s_intg start_of_ped=0;
  s_intg found_per;
  s_intg found_ped=0;
/*s_intg good_format;*/
  s_intg count;             /* Number of people assigned loops.  */
  s_byte done;

  count = 0;

  Rprintf("\n\n\tEnter identifiers for ");
  Rprintf("each pedigree and person...\n");
  Rprintf("\tenter pedigree 0 when finished.\n");

  done = FALSE;
  while(!done) {

    read_pedigree(pedigree_s);

    if (pedigree_s[0] == '0' && pedigree_s[1] == '\0') done = TRUE;
     else {

       found_ped = FALSE;
       i = 1;
	 while ((!found_ped) && (i<=totperson)) {

             if (( strcmp(pedigree_s,person[i]->oldped_s) == 0)) {
	       found_ped = TRUE;
	       start_of_ped = i;
	       pedigree = person[i]->ped;
	     }
           
	   i++;
         }

       if ((i>= totperson) && (!found_ped))
	 Rprintf("\tPedigree not found...\n");
       }

       if((!done) && (found_ped)) {
        found_per = FALSE;
        i = start_of_ped;
        while(!found_per) {

         read_person(person_s);

	 while( (i<=totperson) &&
	       (pedigree == person[i]->ped) &&
	       (!found_per)){
		   if (!strcmp(person[i]->oldid_s,person_s)) {
             loops[count++] = i;
             /* Any loops with indices above i should be incremented by one,
                since add_loop shifts all the indices above i */
             for (ii = 0; ii < count; ii++) 
              if (loops[ii] > i) loops[ii] += 1;
		     found_per = TRUE;
             add_loop(start_of_ped,i);
		   }

		 if (!found_per) i++;
	       } /* end of while */

	 /* The whole pedigree has been searched.  If no match  */
	 /* was found for the identifier then set i back to the */
	 /* begining of the pedigree and re-prompt the user for  */
	 /* for a new id.                                       */

	 if (!found_per) {
	   Rprintf("\tPerson not found...\n");
	   i = start_of_ped;
	 }
       }  /* end of while(!found_per) */
       } /* end of if(!done ) */
  }   /* end of while(!done) */

  save_loops(count);
}
/*****************************************************************************/
/*                                                                           */
/*                          get_loops                                        */
/*                                                                           */
/*****************************************************************************/

#ifdef executable
void get_loops()
#else
void get_loops(int *withloop, char **loopfile)
#endif
{
#ifdef executable
  s_byte response;
  Rprintf("\n");
  Rprintf("Does your pedigree file contain any loops?    (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) {
    Rprintf("\nDo you have a file of loop assignments?       (y/n) -> ");
    fscanf(stdin,"%1s",&response);
    if ((response == 'y') || (response == 'Y')) file_loops();
    else some_loops();
  }
#else
  if(*withloop) file_loops(loopfile);
  else some_loops();
#endif
}
/****************************************************************************/
/*                                                                          */
/*			count_generations				    */
/*                                                                          */
/****************************************************************************/

s_intg count_generations(person_x)
  s_intg person_x;
{
  struct ind *decendent;
  s_intg count;

  count = 0;
  if ( person[person_x]->foff != NULL ) {
    count++;
    decendent = person[person_x]->foff;
    while (decendent->foff != NULL) {
      count++;
      decendent = decendent->foff;
    }
  }
  return(count);
}

/****************************************************************************/
/*                                                                          */
/*			   clear_proband   				    */
/*                                                                          */
/* Given the index of some person, this function clears all proband fields  */
/* in the same pedigree.                                                    */
/*                                                                          */
/****************************************************************************/

void clear_proband(person_index)
     s_intg person_index;
{

  s_intg pedigree;
  s_intg i;
  s_intg found_ped;
  s_intg end_of_ped;

  pedigree = person[person_index]->ped;

  found_ped = FALSE;
  i=1;
  while((!found_ped) && (i<=totperson))
   if (pedigree == person[i++]->ped) found_ped = TRUE;

  if (!cleared[pedigree]){
    i--;
    end_of_ped = FALSE;
    while ((!end_of_ped) && (i<=totperson)) {
      if (pedigree == person[i]->ped) {
	if (person[i]->proband == 1)  /* don't clear a loop person, only clear probands */
	  person[i]->proband = 0;
	else if (person[i]->proband == -1)
	  person[i]->proband = 2;  /* Selecting a proband in this pedigree */
	                           /* => set person in first loop back into loop */
      }
      else end_of_ped = TRUE;
      i++;
    }
    cleared[pedigree] = TRUE;
  }
}

/*****************************************************************************/
/*                                                                           */
/*                          save_probands                                    */
/*                                                                           */
/*****************************************************************************/

void save_probands(count)
    s_intg count;
{
  s_intg i;
  s_byte response;
  s_byte proband_file[max_filespec];
  FILE   *prof;


  Rprintf("\n\nDo you want these selections saved ");
  Rprintf("for later use?  (y/n) -> ");
  fscanf(stdin,"%1s",&response);

  if ((response == 'y') || (response == 'Y')) {
    proband_file[0] = '\0';
    Rprintf("\nEnter filename -> ");
    while ( proband_file[0] == '\0' ) {
      fgets(proband_file,max_filespec,stdin);
    }
    if ( (prof = fopen(proband_file,"w")) == NULL) {
      REprintf("\nERROR: Cannot open file %s\n",proband_file);
    }
    else {
      for(i=0; i<count; i++) {
         fprintf(prof,"%s ",person[ probands[i] ]->oldped_s);
         fprintf(prof,"%s\n",person[ probands[i] ]->oldid_s);
      }
      fclose(prof);
    }
  }
 }

/****************************************************************************/
/*                                                                          */
/*			   auto_probands   				    */
/*                                                                          */
/****************************************************************************/

void auto_probands()
{
  s_intg i,istart;
  s_intg ped_num;		        /* current pedigree */
  s_intg found=0;               /* index of best choice so far */
  s_intg max_level;             /* max num generations found so far */
  s_intg trys;

/* For each male with unknown parents, count the number */
/* of generations which are below him.                  */

  for(i=1; i<=totperson; i++)
   if ((person[i]->paid == 0) &&
       (person[i]->maid == 0) &&
       (person[i]->sex  == 1)) {
       person[i]->generations = count_generations(i);
       }

/* For each pedigree find the first male that has unknown */
/* parents and more generations below him than any other  */
/* male with unknown parents.                             */
i = 1;
istart = 1;
trys = 0;
try_again:
 i = istart;
trys++;
if (trys > 20) {
 REprintf("\nERROR: auto_probands is unable to find in 20 tries a male proband");
 REprintf("\nwho has no parents in the pedigree and");
 REprintf("\nwho is either not in a loop or is in the first loop.");
 REprintf("\nThe problem is in pedigree %s.",person[i]->oldped_s);
 REprintf("\nChange the order in which you choose the loops.\n");
 error("%d",1);
}
while(i<=totperson) {
  ped_num = person[i]->ped;
  max_level = 0;
  istart = i; /* Starting number of this pedigree */  
  while((i<=totperson) && (ped_num == person[i]->ped)) {
    if ( (person[i]->paid == 0) &&
         (person[i]->maid == 0) &&
         (person[i]->sex  == 1) &&
         (person[i]->generations > max_level)) {
      found = i;
      max_level = person[i]->generations;
    }
    i++;
  } /* end of while */
  if (person[found]->proband > 2) {
    /* The person we found is a loopperson, but is not in the first loop */
    /* so we can't use this person as the proband => set generations to -1 */
    /* and try again in this pedigree. */
     person[found]->generations = -1;
     goto try_again;
   } else if (person[found]->proband == 2)
        person[found]->proband = -1; /* Proband chosen is in first loop */
   else person[found]->proband = 1;
 }
}
/*****************************************************************************/
/*                                                                           */
/*                             file_probands                                 */
/*                                                                           */
/*****************************************************************************/
#ifdef executable
void file_probands()
#else
void file_probands(char **probandfile)
#endif
{
/*s_byte response;*/
/*s_intg pedigree;*/
/*s_intg pedigree_i;*/
  s_byte pedigree_s[maxname];
/*s_intg person_i;*/
  s_byte person_s[maxname];
  s_intg /*i,*/j/*,k*/;
/*s_intg from_file;*/
/*s_intg start_of_ped;*/
  s_intg found;
/*s_intg good_format;*/
/*s_byte c;*/
  FILE   *prof;
#ifdef executable 
  s_byte proband_file[max_filespec];
  proband_file[0] = '\0';
  Rprintf("\nEnter filename -> ");
  while ( proband_file[0] == '\0' ) {
       fgets(proband_file,max_filespec,stdin);
  }
  if ( (prof = fopen(proband_file,"r")) == NULL) {
      error("\nERROR: Cannot open file %s\n",proband_file);
  }
#else
   if ((prof = fopen(*probandfile,"r")) == NULL) {
       error("\nERROR: Cannot open file %s\n",*probandfile);
   }
#endif
    auto_probands();  /* makes sure a proband is set for each pedigree */
                      /* even if the input file does not set it.       */
    while(!feof(prof)) {

       fscanf(prof,"%s",pedigree_s);
       fscanf(prof,"%s",person_s);

       if (!feof(prof)) {
       found = FALSE;
       j = 1;
       while(( j<=totperson) && (!found)) {

        if ( (strcmp(pedigree_s,person[j]->oldped_s) == 0) &&
             (strcmp(person_s,person[j]->oldid_s ) == 0)) {
            clear_proband(j);
            if (person[j]->proband > 2) {
   REprintf("\nERROR: If a loopperson is also the proband, that person \n");
   REprintf("       must be in the first loop (#2). \n");
   REprintf("Proband %s in pedigree %s is in loop %d \n"
   ,person[j]->oldid_s,person[j]->oldped_s,person[j]->proband);
   error("%d",1);
   } 
   else person[j]->proband = 1;
   found = TRUE;
        }
        j++;
	if (( j>totperson) && (!found)) {
             error("\nERROR: Ped: %s Per: %s - Not found, check proband file.\n",
		  pedigree_s,
		  person_s);
	}
       }
      }
     }
     fclose(prof);
}



/*****************************************************************************/
/*                                                                           */
/*                              all_probands                                 */
/*                                                                           */
/*****************************************************************************/

void all_probands()
{
/*s_byte response;*/
  s_intg pedigree;
  s_byte person_s[maxname];
/*s_intg person_i;*/
  s_intg i/*,j,k*/;
  s_intg start_of_ped;
  s_intg found;
/*s_intg good_format;*/
  s_intg count;               /* number of people assigned probands.  */

  count = 0;

      Rprintf("\n\nEnter the identifier of the ");
      Rprintf("person who is to be the proband for...\n\n");
      pedigree = 0;
      for(i=1; i<=totperson; i++) {
        if (pedigree != person[i]->ped) {
          start_of_ped = i;
          pedigree = person[i]->ped;

          Rprintf("\n\n\tPedigree   -> ");
          Rprintf("%s\n",person[i]->oldped_s);

          read_person(person_s);

          found = FALSE;
          while( (i<=totperson) &&
                 (pedigree == person[i]->ped) &&
                 (!found)){

             if (!strcmp(person[i]->oldid_s,person_s)) {
              if (person[i]->proband > 2) {
   REprintf("\nERROR: If a loopperson is also the proband, that person \n");
   REprintf("       must be in the first loop (#2). \n");
   REprintf("Proband %s in pedigree %s is in loop %d \n"
   ,person[i]->oldid_s,person[i]->oldped_s,person[i]->proband);
   error("%d",1);
   } else person[i]->proband = 1;
               probands[count++] = i;
               found = TRUE;
             }

             if (!found) i++;
           } /* end of while */

             /* The whole pedigree has been searched.  If no match  */
             /* was found for the identifier then set i back to the */
             /* begining of the pedigree and re-prompt the user for  */
             /* for a new id.                                       */

             if (!found) {
               Rprintf("\tPerson not found...\n");
               i = start_of_ped -1;  /* -1 because the for loop will +1 */
               pedigree = 0;
             }
        }
      }

  save_probands(count);
}

/*****************************************************************************/
/*                                                                           */
/*                              some_probands                                */
/*                                                                           */
/*****************************************************************************/

void some_probands()
{
/*s_byte response;*/
/*s_intg person_i;*/
  s_byte person_s[maxname];
/*s_intg pedigree_i;*/
  s_byte pedigree_s[maxname];
  s_intg pedigree=0;
  s_intg i/*,j,k*/;
  s_intg start_of_ped=0;
  s_intg found_per;
  s_intg found_ped=0;
/*s_intg good_format;*/
  s_intg count;             /* number of people assigned probands.  */
  s_byte done;

  auto_probands();      /* Set all probands before the user selects a few. */
  count = 0;            /* Count number of successful entries. */

  Rprintf("\n\n\tEnter identifiers for ");
  Rprintf("each pedigree and person...\n");
  Rprintf("\tenter pedigree 0 when finished.\n");

  done = FALSE;
  while(!done) {

    read_pedigree(pedigree_s);

    if (pedigree_s[0] == '0' && pedigree_s[1] == '\0') done = TRUE;
     else {

       found_ped = FALSE;
       i = 1;
	 while ((!found_ped) && (i<=totperson)) {
		if (( strcmp(pedigree_s,person[i]->oldped_s) == 0)) {
	       found_ped = TRUE;
	       start_of_ped = i;
	       pedigree = person[i]->ped;
	     }
           
	   i++;
         }

       if ((i>= totperson) && (!found_ped))
	 Rprintf("\tPedigree not found...\n");
       }

       if((!done) && (found_ped)) {
        found_per = FALSE;
        i = start_of_ped;
        while(!found_per) {

          read_person(person_s);

	 while( (i<=totperson) &&
	       (pedigree == person[i]->ped) &&
	       (!found_per)){
		   if (!strcmp(person[i]->oldid_s,person_s)) {
		     clear_proband(i);
              if (person[i]->proband > 2) {
   REprintf("\nERROR: If a loopperson is also the proband, that person \n");
   REprintf("       must be in the first loop (#2). \n");
   REprintf("Proband %s in pedigree %s is in loop %d \n"
   ,person[i]->oldid_s,person[i]->oldped_s,person[i]->proband);
   error("%d",1);
   } else person[i]->proband = 1;
             probands[count++] = i;
		     found_per = TRUE;
		   }

		 if (!found_per) i++;
	       } /* end of while */

	 /* The whole pedigree has been searched.  If no match  */
	 /* was found for the identifier then set i back to the */
	 /* begining of the pedigree and re-prompt the user for  */
	 /* for a new id.                                       */

	 if (!found_per) {
	   Rprintf("\tPerson not found...\n");
	   i = start_of_ped;
	 }
       }  /* end of while(!found_per) */
       } /* end of if(!done ) */
  }   /* end of while(!done) */

  save_probands(count);
}

/*****************************************************************************/
/*                                                                           */
/*                          get_probands                                     */
/*                                                                           */
/*****************************************************************************/
#ifdef executable
void get_probands()
#else
void get_probands(int *auto_proband,char **probandfile)
#endif
{
  s_byte response;
  response='y';
#ifdef executable
  Rprintf("\n");
  Rprintf("Do you want probands selected automatically?   (y/n) -> ");
  fscanf(stdin,"%1s",&response);
#else
  if (*auto_proband) response='y';
#endif
  if ((response == 'y') || (response == 'Y')) auto_probands();
  else {
#ifdef executable
  Rprintf("\nDo you have a file of proband assignments?    (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) file_probands();
#else
  response='y';
  if ((response == 'y') || (response == 'Y')) file_probands(probandfile);
#endif
  else {
 
  Rprintf("\nDo you want to select all probands?           (y/n) -> ");
  fscanf(stdin,"%1s",&response);
  if ((response == 'y') || (response == 'Y')) all_probands();
  else some_probands();
}
}
}

/****************************************************************************/
/*                                                                          */
/*                           writeped                                       */
/*                                                                          */
/****************************************************************************/

void writeped()
{
  s_intg i;
  s_byte *ped_format;
  s_byte *ind_format;
  s_byte *format_1  = "%1d";
  s_byte *format_2  = "%2d";
  s_byte *format_3  = "%3d";
  s_byte *format_4  = "%4d";
  s_byte *format_5  = "%5d";
  s_byte *format_6  = "%6d";
  s_byte *format_7  = "%7d";
  s_byte *format_8  = "%8d";
  s_byte *format_9  = "%9d";
  s_byte *format_10 = "%10d";

  /* setup pedigree printing format */

  if (biggest_p_id >= 10e10)
    ped_format = format_10;
  else
  if (biggest_p_id >= 10e09)
    ped_format = format_9;
  else
 if (biggest_p_id >= 10e08)
    ped_format = format_8;
  else
  if (biggest_p_id >= 10e07)
    ped_format = format_7;
  else;
  if (biggest_p_id >= 10e06)
    ped_format = format_6;
  else
  if (biggest_p_id >= 10e05)
    ped_format = format_5;
  else
  if (biggest_p_id >= 10e04)
    ped_format = format_4;
  else
  if (biggest_p_id >= 10e03)
    ped_format = format_3;
  else
  if (biggest_p_id >= 10e02)
    ped_format = format_2;
  else
    ped_format = format_1;

  /* setup individual printing format */

  if (biggest_i_id > 9999)
    ind_format = format_6;
  else
  if (biggest_i_id > 999)
    ind_format = format_5;
  else
  if (biggest_i_id > 99)
    ind_format = format_4;
  else
  if (biggest_i_id > 9)
    ind_format = format_3;
  else
    ind_format = format_2;


  for( i=1; i<=totperson; i++) {
    if (ped_integers)
    fprintf(pedout,"%s",  person[i]->oldped_s);     /* Use original ids */
    else fprintf(pedout,ped_format, person[i]->ped); /* This renumbers pedigrees */
    
    fprintf(pedout,ind_format,  person[i]->id);

    if (person[i]->pa != NULL)
      fprintf(pedout,ind_format,  person[i]->pa->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->ma != NULL)
      fprintf(pedout,ind_format,  person[i]->ma->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->foff != NULL)
      fprintf(pedout,ind_format,  person[i]->foff->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->nextpa != NULL)
      fprintf(pedout,ind_format,  person[i]->nextpa->id);
    else
      fprintf(pedout,ind_format,  0);

    if (person[i]->nextma != NULL)
      fprintf(pedout,ind_format,  person[i]->nextma->id);
    else
      fprintf(pedout,ind_format,  0);

    fprintf(pedout,"%2d", person[i]->sex);
    /* If we never selected a proband in this pedigree, then the       */
    /* person in the first loop selected by autoprob remains the proband */
    if (person[i]->proband == -1) person[i]->proband = 1;
    fprintf(pedout,"%2d", person[i]->proband);
    fprintf(pedout,"%s", person[i]->phen->phen_chars);
    fprintf(pedout,"  Ped: %s",person[i]->oldped_s);
    fprintf(pedout,"  Per: %s\n",person[i]->oldid_s);
   }
}

/*****************************************************************************/
/*                                                                           */
/*                             strcmp_i                                      */
/*                                                                           */
/* A case insensitive string compare. Returns 0 if equal, 1 if not.          */
/*                                                                           */
/*****************************************************************************/

s_intg strcmp_i(s1,s2)
     register u_byte *s1;
     register u_byte *s2;
{
  while(((*s1 >= 'a' && *s1 <= 'z')?
	 (*s1 & 0xdf):
	 (*s1))==
	((*s2 >= 'a' && *s2 <= 'z')?
	 (*s2++ & 0xdf):
	 (*s2++)))
    if (*s1++ == '\0')
      return(0);
  return(1);
}

/*****************************************************************************/
/*                                                                           */
/*                               check_ids                                   */
/*                                                                           */
/*****************************************************************************/

void check_ids()
{
  s_intg i,j;
  
  /* Scan through the pedigree ids until the first non-digit is found; */
  /*  set ped_integers to FALSE if any non-digit is found.             */ 
  ped_integers = TRUE;
  i = 1;
  while ((ped_integers) && (i<=totperson)) {
    for (j=0; person[i]->oldped_s[j] != '\0'; j++) 
    if (!isdigit(person[i]->oldped_s[j])) {
       ped_integers = FALSE;
       break;
    }
    i++;
  }
}

/*****************************************************************************/
/*                                                                           */
/*                               check_sex                                   */
/*                                                                           */
/*****************************************************************************/

void check_sex()
{
  s_intg i;

  for(i=1; i<=totperson; i++) {

    /* Verify that each person has either 0 or 2 parents. */

    if(((person[i]->pa != NULL) && (person[i]->ma == NULL)) ||
       ((person[i]->pa == NULL) && (person[i]->ma != NULL))){
       REprintf("\nERROR: Ped: %s  Per: %s - Only one parent.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
    }

    /* Verify that father is male. */

   if((person[i]->pa != NULL) && (person[i]->pa->sex != 1)) {
       REprintf("\nERROR: Ped: %s  Per: %s - Sex of father.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
    }

    /* Verify that mother is female. */

   if((person[i]->ma != NULL) && (person[i]->ma->sex != 2)) {
       REprintf("\nERROR: Ped: %s  Per: %s - Sex of mother.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
    }
  }
}

/*****************************************************************************/
/*                                                                           */
/*                           check_no_phen                                   */
/*                                                                           */
/* If a person has no phenotypic data then they were referenced as a parent  */
/* in the data file but did not have a record of their own.                  */
/*                                                                           */
/*****************************************************************************/

void check_no_phen()
{

  s_intg i;

  for(i=1; i<=totperson; i++) {

   if (person[i]->phen == NULL){
       REprintf("\nERROR: Ped: %s  Per: %s - No data.\n",
	      person[i]->oldped_s,person[i]->oldid_s);
       found_error = TRUE;
   }
 }
}

/*****************************************************************************/
/*                                                                           */
/*                           check_no_family                                 */
/*                                                                           */
/* If a person has no parents and no children then report it.                */
/*                                                                           */
/*****************************************************************************/

void check_no_family()
{
  u_intg i;

  /* Mark each person that is referenced as a parent of another person. */

  for( i=1; i<=totperson; i++) {
    if (person[i]->pa != NULL) person[i]->pa->is_parent = TRUE;
    if (person[i]->ma != NULL) person[i]->ma->is_parent = TRUE;
  }

  /* Report all persons who do not have their is_parent  */
  /* flag set and do not have any parents of their own.  */

  for( i=1; i<=totperson; i++ ) {
    if ((!person[i]->is_parent) &&
        (person[i]->pa == NULL) &&
        (person[i]->ma == NULL)) {
	REprintf("\nERROR: Ped: %s  Per: %s - No family.\n",
		person[i]->oldped_s,person[i]->oldid_s);
	found_error = TRUE;
	}
  }
}

#ifdef executable
/*****************************************************************************/
/*                                                                           */
/*                                 main                                      */
/*                                                                           */
/*****************************************************************************/

main(argc,argv)
     s_intg argc;
     u_byte *argv[];
{

  s_intg no_ques = FALSE;
  u_byte response;

  Rprintf("\n           MAKEPED Version %4.2f\n\n",MAKEPED_VERSION);
  Rprintf("Usage: makeped pedigree.file output.file n \n");
  Rprintf("where the 'n' argument tells MAKEPED that \n");
  Rprintf("  1) there are no loops \n");
  Rprintf("  2) all probands should be chosen automatically \n\n");
  Rprintf("This version treats all ids as strings\n");
  Rprintf("but will use original pedigree ids if they ALL are integers\n");
  Rprintf("In addition, this version allows for multiple loops\n");
/*  fprintf(stdout, "Please report all problems and suggestions to \n");
  Rprintf("  Daniel E. Weeks or Jurg Ott \n");
  Rprintf("  Columbia University, Unit 58 \n");
  Rprintf("  722 West 168th Street  \n");
  Rprintf("  New York, NY 10032 \n");
  Rprintf("  Tel (212) 960-2507  Fax (212) 568-2750 \n\n"); */
  
  Rprintf(" Constants in effect \n");
  Rprintf("Maximum number of pedigrees                %d\n",maxped-1);
  Rprintf("Maximum number of individuals              %d\n",maxind-1);
  Rprintf("Maximum characters used in phenotypic data %d\n", maxallchars);
  Rprintf("Maximum number of characters in an id      %d\n\n",maxname);
 
  if (argc > 4) {
    error("\nERROR: Two many command line arguments");
  }

  pifile[0] = '\0';
  pofile[0] = '\0';

  if( argc > 1) {                        /* GET FILESPEC IF ON COMMAND LINE */
    strcpy(pifile, argv[1]);
  }

  else{                                 /* FILES ARE NOT ON COMMAND LINE */
    while ( pifile[0] == '\0' ) {
      Rprintf("Pedigree file -> ");
      fgets(pifile,max_filespec,stdin);
    }
    while ( pofile[0] == '\0') {
      Rprintf("Output file   -> ");
      fgets(pofile,max_filespec,stdin);
    }
  }


  if( argc > 2) {                       /* GET OUTPUT FILE IF ON COM LINE */
    strcpy(pofile, argv[2]);
  }
  else {
    while ( pofile[0] == '\0') {
      Rprintf("Output file   -> ");
      fgets(pofile,max_filespec,stdin);
    }
  }
  /* If the third argument is 'n', then ask no questions -> there  */
  /* are no loops and the probands are to be selected automatically */
  if ( argc > 3) {
    if (*argv[3] == 'N' || *argv[3] == 'n')
      no_ques = TRUE;
  }
  found_error = FALSE;

  if ((pedfile = fopen(pifile, "r")) == NULL){
   error("\nERROR: Cannot open %s\n",pifile);
 }

  if ((pedout = fopen(pofile, "w")) == NULL){
   error("\nERROR: Cannot open %s\n",pofile);
 }
  readped();
  check_ids();
  check_sex();
  check_no_phen();
  check_no_family();
  if(found_error) error("%d",1);
  pointers();
  if (no_ques) auto_probands(); /* Set all probands automatically */
  else {
   get_loops(); 
   get_probands();
  }
  writeped();

  fclose(pedfile);
  fclose(pedout);

}
#else 
void makeped_c(char **pifile, char **pofile, int *autoselect, 
             int *withloop, char **loopfile, 
             int *autoproband, char **probandfile)
{
  Rprintf("\n           MAKEPED Version %4.2f\n\n",MAKEPED_VERSION);
  Rprintf(" Constants in effect \n");
  Rprintf("Maximum number of pedigrees                %d\n",maxped-1);
  Rprintf("Maximum number of individuals              %d\n",maxind-1);
  Rprintf("Maximum characters used in phenotypic data %d\n", maxallchars);
  Rprintf("Maximum number of characters in an id      %d\n\n",maxname);
 
  found_error = FALSE;

  if ((pedfile = fopen(*pifile, "r")) == NULL){
   error("\nERROR: Cannot open %s\n",*pifile);
 }

  if ((pedout = fopen(*pofile, "w")) == NULL){
   error("\nERROR: Cannot open %s\n",*pofile);
 }
  readped();
  check_ids();
  check_sex();
  check_no_phen();
  check_no_family();
  if(found_error) error("%d",1);
  pointers();
  if (*autoselect) auto_probands(); /* Set all probands automatically */
  else {
   get_loops(withloop,loopfile);
   get_probands(autoproband,probandfile);
  }
  writeped();

  fclose(pedfile);
  fclose(pedout);

}

#endif

/* Adapted from makeped.c on 18-10-2003 */
/* file_loops Change u_byte to s_byte pedigree_s on 16-8-2004 */
/* Change gets to fgets in relation to console input on 2-12-2013 */
