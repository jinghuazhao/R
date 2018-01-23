/*21/4/2004 adapt for GAW14 data problem 1 (COGA data)

  Usage:
  1. unzip GAW14 data files, then the data will be in Microsat directory
  2. compile gaw14.c
  3. issue command gaw14

  Changes:
  . new file names
  . more pedigrees, longer directory names, wider fields, more information
  . sex in pheno.dat is now character
  . no pedigree id in .dat file and the xxx/ xxx format changed to xxx/xxx
  . pedigrees 27 and 106 now have loops, confirmed by kinship() and loops
  . follow LINKAGE to use filenames .dat, .pre, .ped and .lop
  . comment on these changes!
*/

/*Utility program for GAW11 problem 1*/
/*draft on 23/6/1998 JH Zhao*/
/*add disease model on 2/7/1998 (also became restricted to Unix system only)*/
/*add member clips for genehunter on 20/7/1998*/
/*to separate gaw14.h on 16/5/2004*/

#include "gaw14.h"

main(int argc,char **argv)
{
int i,j,l,chromsom,totloci;
/**/
char parfile[25],pedfile[25],frqfile[25],mapfile[25],datfile[25],makfile[25],ghfile[25],ghout[25];
ind p;
LOC q;

person=(ind *)calloc(MAX_OBS,sizeof(ind));
getphen();
totloci=0;
for(chromsom=1;chromsom<24;++chromsom)
{
  if(chromsom<10)
  {
    sprintf(parfile,"chr0%-d.dat",chromsom);
    sprintf(pedfile,"chr0%-d.pre",chromsom);
    sprintf(frqfile,"Microsat/ms0%-d.frq",chromsom);/**/
    sprintf(mapfile,"Microsat/ms0%-d.map",chromsom);
    sprintf(datfile,"Microsat/ms0%-d.dat",chromsom);
    sprintf(makfile,"chr0%-d.mak",chromsom);
    sprintf(ghfile,"chr0%-d.in",chromsom);
    sprintf(ghout,"chr0%-d.out",chromsom);
  }
  else
  {
    sprintf(parfile,"chr%-d.dat",chromsom);
    sprintf(pedfile,"chr%-d.pre",chromsom);
    sprintf(frqfile,"Microsat/ms%-d.frq",chromsom);
    sprintf(mapfile,"Microsat/ms%-d.map",chromsom);
    sprintf(datfile,"Microsat/ms%-d.dat",chromsom);
    sprintf(makfile,"chr%-d.mak",chromsom);
    sprintf(ghfile,"chr%-d.in",chromsom);
    sprintf(ghout,"chr%-d.out",chromsom);
  }
  par=fopen(parfile,"w+");
  if(par==NULL) perror("Error opening LINKAGE parameter file !");
  getfrq(frqfile);totloci+=nloci;
  getmap(mapfile);
  fclose(par);

  ped=fopen(pedfile,"w+");
  if(ped==NULL) perror("Error opening LINKAGE data file !");
  genotype=(LOC *)calloc(MAX_OBS,sizeof(LOC));
  getdat(datfile);
  a2a();
  l=0;
  for(i=0;i<sample_size;++i)
  {
    p=person[i];
/*  if(p.id==clips[l]) {l++;continue;}*/ /**/
    fprintf(ped,"%6d%9d%9d%9d%2s%2d",p.pid,p.id,p.fid,p.mid,(strcmp(p.sex,"F")?"1":"2"),p.aff);
    q=genotype[i];
    for(j=0;j<nloci;++j) fprintf(ped,"  %2d %2d",q.allele[j][0],q.allele[j][1]);
    fprintf(ped,"\n");
  }
  free(genotype);
  fclose(ped);

  makeped=fopen(makfile,"w+");
  if(makeped==NULL) perror("Error opening MAKEPED command file !");
  if(chromsom<10)
  {
    fprintf(makeped,"chr0%-d.pre\n",chromsom);
    fprintf(makeped,"chr0%-d.ped\n",chromsom);
  }
   else
  {
    fprintf(makeped,"chr%-d.pre\n",chromsom);
    fprintf(makeped,"chr%-d.ped\n",chromsom);
  }
  fprintf(makeped,"y\n"); /*loops==yes*/
  fprintf(makeped,"y\ngaw14.lop\n"); /*a file==yes*/
  fprintf(makeped,"y\n"); /*proband*/
  fclose(makeped);

  gh=fopen(ghfile,"w+");
  if(gh==NULL) perror("Error opening Genehunter command file !");
  fprintf(gh,"system rm -f %s \n",ghout);
  fprintf(gh,"photo %s \n",ghout);
  fprintf(gh,"time \n");
  fprintf(gh,"load %s \n",parfile);
  fprintf(gh,"use ");
  for(i=1;i<nloci;++i) fprintf(gh,"%d %.1f ",i,map[i]);fprintf(gh,"%d \n",nloci);
  fprintf(gh,"skip large off \n");
  fprintf(gh,"single on \n");
  fprintf(gh,"scan %s \n",pedfile);
  fprintf(gh,"total stat het \n");
  fprintf(gh,"time \n");
  fprintf(gh,"quit \n");
  fclose(gh);

}
free(person);
lop=fopen("gaw14.lop","w");
if(!lop)
{
  fprintf(stderr,"I cannot open output for loop file gaw14.lop\n");
  exit(1);
}
fprintf(lop,"10051 10000529\n");
fprintf(lop,"10065 10001161\n");
fclose(lop);
#ifdef DEBUG
  printf("%5d individuals in %3d families and %d loci \n",sample_size,families,totloci);
  printf("Family ID Size\n");
  for(i=0;i<families;++i) printf("%4d %4d \n",i+1,fsize[i]);
#endif
return 0;
}

int getfrq(char *locfile)
/*read marker allele frequencies*/
{
FILE *fp;
char line[LINELEN],name[10];
mark m;
int i,j,size;
float freq,s;

fp=fopen(locfile,"r");
if(!fp) fprintf(stderr,"cannot open file %s\n",locfile);
j=0;
while(fgets(line,LINELEN,fp)&&sscanf(line,"%s %d",&m.name,&m.n)==2)
{
   s=0;
   for(i=0;i<m.n;++i)
   {
      fgets(line,LINELEN,fp);
      sscanf(line,"%d %f",&m.size[i],&m.freq[i]);
      s+=m.freq[i];
   }
   locus[j]=m;++j;
   if(fabs(s-1.0)>=0.0001) printf("%10s has allele frequencies summed up to %f \n",m.name,s);
}
fclose(fp);
nloci=j;

fprintf(par,"%d  0 0 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM\n",nloci+1);
fprintf(par," 0 0.0 0.0 0  << MUT LOCUS, MUT RATE, HAPLOTYPE FREQUENCIES (IF 1)\n");
for(i=0;i<nloci+1;++i) fprintf(par,"%3d ",i+1);fprintf(par,"\n");
fprintf(par,"1   2  << AFFECTION, NO. OF ALLELES\n");
fprintf(par," %.6f %.6f   << GENE FREQUENCIES\n",1-0.014,0.014);
fprintf(par," 1 << NO. OF LIABILITY CLASSES\n");
fprintf(par," %.4f %.4f %.4f << PENETRANCES\n",0.04,0.4,0.4);
for(j=0;j<nloci;++j)
{
   m=locus[j];
   fprintf(par,"3 %2d    # %-10s %15d\n",m.n,m.name,j+1);
   for(i=0;i<m.n;++i) fprintf(par,"%.4f ",m.freq[i]);
   fprintf(par,"\n");
}
fprintf(par," 0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)\n");
return 0;
}

int getmap(char *mapfile)
/*read maps*/
{
FILE *fp;
char line[100],rest[100],name[10];
float di,d[MAX_LOC];
int i,j;

fp=fopen(mapfile,"r");
i=0;
while(fgets(line,100,fp)&&sscanf(line,"%s %f %[^\n]",&name,&di,rest)>1) d[i++]=di;
fclose(fp);
for(i=1;i<nloci;++i) map[i]=d[i]-d[i-1];
fprintf(par,"0.5 ");
for(j=1;j<i;++j) fprintf(par,"%.6f ",haldane(map[j]/100));fprintf(par,"\n");
fprintf(par," 1 0.10000 0.45000 << REC VARIED, INCREMENT, FINISHING VALUE\n");
return 0;
}

int getdat(char *datfile)
/*read genotype data*/
{
FILE *fp;
char line[LINELEN],rest[LINELEN],a1[5],a2[5],*ptr;
int i,j;
LOC p;

fp=fopen(datfile,"r");
i=0;
while(fgets(line,LINELEN,fp)&&sscanf(line,"%d%[^\n]",&p.id,rest)>1)
{
   strcpy(line,rest);
   for(j=0;j<nloci;++j,strcpy(line,rest),*rest='\0')
   {
      p.allele[j][0]=p.allele[j][1]=0;
      ptr=line;
      strncpy(a1,ptr,4);sscanf(a1,"%d",&p.allele[j][0]);/*1-4*/
      strncpy(a2,ptr+5,3);sscanf(a2,"%d",&p.allele[j][1]);/*5th*/
      strncpy(rest,ptr+8,strlen(line)-8);/*jump width=8*/
   }
   genotype[i]=p;
   ++i;
}
fclose(fp);
return 0;
}

int a2a()
/*assign allele numbers, assuming increase by allele size*/
{
int i,j,k,l;
LOC p;

for(i=0;i<sample_size;++i)
{
  p=genotype[i];
  for(j=0;j<nloci;++j)
  {
     for(k=0;k<locus[j].n;++k) for(l=0;l<2;++l)
     {
        if(p.allele[j][l]==locus[j].size[k])
        {
          p.allele[j][l]=k+1; continue;
        }
     }
}
genotype[i]=p;
}
return 0;
}

int getphen()
/*read phenotype data*/
{
FILE *fp;
char line[LINELEN],rest[LINELEN];
int i,j,k,l,l1,l2,id,age,race,aff;
ind p;

fp=fopen("pheno.dat","r");
if(fp==NULL) perror("Error opening pheno.dat !");/**/
i=0;
while(fgets(line,LINELEN,fp))
{
   sscanf(line,"%6d%9d%9d%9d%[^\n]",&p.pid,&p.id,&p.fid,&p.mid,rest);/**/
   strcpy(line,rest);
   sscanf(line,"%2s%*d%*d%*d%3d%2d%d%[^\n]",&p.sex,&age,&race,&aff,rest);/**/
   strcpy(line,rest);
   p.aff=0;
   if(age>25&&aff==1) p.aff=1;
   if(aff==5) p.aff=2;
   person[i]=p;
   ++i;
}
fclose(fp);
sample_size=i;
id=0;k=1;i=0;
l1=person[id].pid;
for(j=0;j<sample_size;++j){
   id++;
   l2=person[id].pid;
   if(l1!=l2&&l1!=0)
   {
     fsize[i++]=k;
     l1=l2;k=1;
   }
   else ++k;
}
families=i;
return 0;
}
