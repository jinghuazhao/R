/* Header file for C programs for nuclear families */

typedef struct offspring {/* List element for siblings */
  int id;                 /* Identifier */
  int affected;           /* 0 = U/k, 1=No, 2=Yes */
  int *markers;           /* 2xm matrix of markers */
  int *ivec_f, *ivec_m;   /* Inheritance vectors */
  int f_tr,  m_tr;        /* Which paternal/maternal haplotype transmitted */
  int sex;                /* sex (optional) */
  double qt;              /* quantitative trait phenotype (optional) */
  struct offspring *next; /* pointer to next sibling */
} Offspring;

typedef struct family {   /* List element for a single nuclear family */
  int check;             /* 3=complete, 2/1 = mother/father only, 0=neither */
  int pedigree, father_id, mother_id; /* Identification data */
  int *father, *mother;   /* Parental haplotypes as 2xm matrices */
  int *phase_f, *phase_m; /* Haplotype phase vectors*/
  Offspring *children;    /* list of children */
  struct family *next;    /* pointer to next family */
} Family;

/* S and R callable functions */

void hap_transmit(int *n, int *ped, int *id, int *father, int *mother,
		  int *sex, int *aff, int *if_qt, double *qt, 
		  int *m, int *markers, 
		  int *multiple_cases, int *impute_using_affected,
		  char **ofname); /* Write haplotypes to disk */

void hap_read(int *n, int *ped, int *id, int *father, int *mother,
	      int *if_qt, double *qt, 
	      int *m, int *f_tr, int *f_un, int *m_tr, int *m_un, 
	      char **ifname); /* Read haplotype file back from disk */

/* Functions */

Family *new_family(int m);   /* Create a family object */
Family *copy_family(Family *f, int m);  /* Copy family object */
void del_family(Family *f);  /* Delete a family object */
int count_families(Family *f); /* Count th efamilies in a list */
int count_offspring(Family *first, int affected_only); /* Count children */
Offspring *new_child(int m);   /* Create an offspring object */
Offspring *copy_child(Offspring *child, int m);  /* Copy offspring object */ 
void del_child(Offspring *child); /* Delete an offspring object */
void show_family(Family *f); /* Summary of family to file */
void print_family(Family *f, int m, FILE *stream); /* Print data */
void warn(char *message, Family *f); /* Warning to stderr */
int inheritance(Family *f, int m); /* Compute inheritance vectors */
int haplotype(Family *f, int m, int resolve_homozygous); 
                                  /* Resolve haplotypes and transmission */
int impute_parent(Family *f, int m, int use_affected);
/* Repeat a nuclear family with one affected per copy */
Family *expand_family(Family *f, int m);
/* Treat second+ affected offspring as unknown status */
void use_only_first(Family *f); 
/* Construct list of nuclear families  from pedfile */
Family *nuclear(int n, int *ped, int *mem, int *father, int *mother,
		int *sex, int *aff_status, double *qt, 
		int m, int *markers);
/* Write a list of transmitted and untransmitted htypes for affected children*/
int hap_write(Family *first, int m, int if_qt, FILE *stream);
