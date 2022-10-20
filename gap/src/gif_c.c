#include <R.h>
#include <stdio.h>
#include <string.h>
#ifdef UNIX
#include <sys/time.h>
#include <sys/resource.h>
#else
#include <time.h>
#include <stdlib.h>
#endif
/*===============================*/
/* Structures for graph handling */
/*===============================*/

typedef struct vertex
{
  int             proband, id;
  struct edge     *up, *down;
  struct vertex   *left, *right;
} vertex ;

typedef struct edge
{
  int             free;
  struct vertex   *top, *bottom;
  struct edge     *next_up, *next_down;
} edge ;

typedef struct vertex_list
{
  struct vertex           *c;
  struct vertex_list      *n;
} vertex_list;

/*================================*/
/* Graph handling routines needed */
/*================================*/

void            no_probands(void), make_edge(vertex *, vertex *);
int             connected(vertex *, vertex *), new_proband(vertex *);
vertex          *find_vertex(int);
vertex_list     *proband_list;

void gif_c(int *data, int *famsize, int *gifset, int *giflen, double *gifval)
{
  int     id, i, j, k, n_prob;
  double  total_kinship(void);
  vertex  *top, *bot;

  top=bot=NULL;
  for (id=0; id<*famsize; id++)
  {
    i=data[id*3];
    j=data[id*3+1];
    k=data[id*3+2];
    if (i > 0) bot = find_vertex(i);
    if (j > 0)
    {
      top = find_vertex(j);
      if (!connected(bot,top)) make_edge(bot,top);
    }
    if (k > 0)
    {
      top = find_vertex(k);
      if (!connected(bot,top)) make_edge(bot,top);
    }
  }
  no_probands();
  n_prob=0;
  for (id=0;id<*giflen;id++)
  {
    i=gifset[id];
    if (i > 0)
    {
      bot=find_vertex(i);
      if (new_proband(bot)) n_prob +=1;
    }
  }
  *gifval=100000.0*total_kinship()/n_prob/(n_prob-1)*2.0;

}

/*==========================*/
/* Kinship finding routines */
/*==========================*/

static int max_path_length=200, path_length[200];
void path_find(vertex *v, int length, int going_up)
{
  edge    *e;
  int     new_length = length+1;
  if (v->proband) path_length[new_length] += 1;
  if (going_up)
  {
    for (e = v->up; e; e = e->next_up)
    {
      e->free = 0;
      path_find(e->top,new_length,1);
      e->free = 1;
    }
  }

  for (e = v->down; e; e = e->next_down) if (e->free)
  {
    e->free = 0;
    path_find(e->bottom,new_length,0);
    e->free = 1;
  }
}

double total_kinship(void)
{
  int             i;
  double          half, kin;
  vertex_list     *vl;
  for (i=0; i<max_path_length; i++) path_length[i] = 0;
  for (vl=proband_list; vl; vl = vl->n)
  {
    vl->c->proband = 0;
    path_find(vl->c,0,1);
  }
  for (i=2, half=0.5, kin=0.0; i<max_path_length; i++)
  {
    half *= 0.5;
    kin += path_length[i] * half;
  }
  return(kin);
}

/*=========================*/
/* Graph handling programs */
/*=========================*/

static vertex   *binary_tree = NULL;
static int      n_proband=0;

vertex  *new_vertex(void)
{
  vertex  *v;
  if (!(v = (vertex *) calloc(1,sizeof(vertex))))
  {
    error("\nnew_vertex: cannot allocate vertex");
  }
  return(v);
}

vertex *find_vertex(int i)
{
  vertex  **vv;

  for (vv = &(binary_tree); *vv && (*vv)->id != i; )
  {
    vv = (i < (*vv)->id ? &((*vv)->right) : &((*vv)->left) );
  }
  if (!(*vv))
  {
    *vv = new_vertex();
    (*vv)->id = i;
  }
  return(*vv);
}

vertex_list *new_vertex_list(void)
{
  vertex_list     *vl;
  if (!(vl = (vertex_list *)calloc(1,sizeof(vertex_list))))
  {
    error("\nnew_vertex_list: cannot alloate vertex_list");
  }
  return(vl);
}

int new_proband(vertex *v)
{
  vertex_list     *vl;

  if (v->proband) return(0);
  v->proband = 1;
  vl = new_vertex_list();
  vl->c = v;
  vl->n = proband_list;
  proband_list = vl;
  n_proband++;
  return(1);
}

edge *new_edge(void)
{
  edge    *e;
  if (!(e = (edge *)calloc(1,sizeof(edge))))
  {
    error("\nnew_edge: cannot allocate edge");
  }
  e->free = 1;
  return(e);
}

void make_edge(vertex *bot, vertex *top)
{
  edge    *e;

  e = new_edge();
  e->top = top;
  e->next_down = top->down;
  top->down = e;
  e->bottom = bot;
  e->next_up = bot->up;
  bot->up = e;
}

int connected(vertex *bot, vertex *top)
{
  edge    *e;

  for (e = bot->up; e && e->top != top; e = e->next_up);
  if (e)  return(1);
  else    return(0);
}

void    free_edge(edge *e)
{
  if (e)
  {
    free_edge(e->next_up);
    free((char *)e);
  }
}

void free_vertex(vertex *v)
{
  if (v)
  {
    free_vertex(v->left);
    free_vertex(v->right);
    free_edge(v->up);
    free((char *)v);
  }
}

void free_vertex_list(vertex_list *vl)
{
  if (vl)
  {
    free_vertex_list(vl->n);
    free((char *)vl);
  }
}

void no_probands(void)
{
  vertex_list *vl;
  if (proband_list)
  {
    for (vl = proband_list; vl; vl=vl->n) vl->c->proband = 0;
    free_vertex_list(proband_list);
    proband_list = NULL;
  }
  n_proband = 0;
}

/*==============*/
/* Time program */
/*==============*/

#ifdef UNIX
double cpu_time(void)
{
  struct rusage   r;
  double          time;

  getrusage(RUSAGE_SELF,&r);
  time = r.ru_utime.tv_sec;
  time += r.ru_utime.tv_usec/1000000.0;
  return(time);
}

#else

double cpu_time (void)
{
  time_t   start, finish;
  double   elapsed_time;

  time( &start );
  time( &finish );

  elapsed_time = difftime( finish, start );
  return elapsed_time;
}
#endif

#ifdef executable

/*===============================*/
/* Parameters for input handling */
/*===============================*/

#define MAXLINE 1000 /* Maximum input line length */
char            *blankline, whereblank;
static char     line_buff[MAXLINE];
int             line_no;

/*==============*/
/* Time program */
/*==============*/

double  cpu_time(void);

/*=========================*/
/* Input handling routines */
/*=========================*/

char *newline()
{
  int i, blank;

  if (!fgets(line_buff,MAXLINE,stdin)) return(NULL);
  while (line_buff[0] == '%')
  {
    if (line_buff[1] == '#') Rprintf("%%# %6.3f (0.000 on PC)\n",cpu_time());
    else if (line_buff[1] =='%') Rprintf("%s",line_buff);
    R_FlushConsole();
    if (!fgets(line_buff,MAXLINE,stdin)) return(NULL);
  }
  for (i=0, blank=1; line_buff[i] != '\n' && blank; i++)
    blank = blank && (line_buff[i] == '\t' || line_buff[i] == ' ');
  if (blank) return(blankline);
  return(line_buff);
}

int got_opt(char *c, char **v, char *o)
/* Returns 1 if string o is one of the c strings in v, and removes it */
{
  int i, j;
  for (i=0; i<*c; i++) if (!strcmp(v[i],(char *)o))
  {
    for (j=i,(*c)--; j<*c; j++) v[j] = v[j+1];
    return(1);
  }
  return(0);
}

void print_profile(int n)
{
  int     i, m;
  m = (n < max_path_length ? n : max_path_length);
  for (i=2; i<m; i++) Rprintf(" %d",path_length[i]);
}

/*==============*/
/*==============*/
/* Main program */
/*==============*/
/*==============*/

int main(int argc, char **argv)
{
  char    *progname = argv[0], *line, *newline();
  int     i, j, k, n_prob, got_opt();
  double  total_kinship(void);
  vertex  *top, *bot;
  blankline = &whereblank;
  line_no=0;
  for (line=newline(); line && line != blankline; line = newline())
  {
    if (sscanf(line,"%d%d%d",&i,&j,&k) != 3)
    {
      error("\n %s(%d): cannot read triplet",progname,line_no);
    }
    if (i > 0) bot = find_vertex(i);
    if (j > 0)
    {
      top = find_vertex(j);
      if (!connected(bot,top)) make_edge(bot,top);
    }
    if (k > 0)
    {
      top = find_vertex(k);
      if (!connected(bot,top)) make_edge(bot,top);
    }
  }

  for ( ; line; line=newline())
  {
    while (line && line == blankline) line=newline();
    for (no_probands(),n_prob=0 ; line && line != blankline; line=newline())
    {
      if (sscanf(line,"%d",&i) != 1)
      {
        error("\n %s(%d): cannot read integer",progname,line_no);
      }
      if (i > 0)
      {
        bot=find_vertex(i);
        if (new_proband(bot)) n_prob +=1;
      }
    }
    Rprintf("%9.3f",100000.0*total_kinship()/n_prob/(n_prob-1)*2.0);
    Rprintf("\n\n");
    R_FlushConsole();
  }
  R_ClearerrConsole();
  return(0);
}

#endif
