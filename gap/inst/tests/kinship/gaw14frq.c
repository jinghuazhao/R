/*16-05-2004 extracted from gaw14.c*/
/*7-6-2004 done for 23 chromosomes*/
/*23-6-2004 rescale allele frequcencies to suppress complaints*/
/*to convert GAW14 frequency format to SOLAR format*/

#include "gaw14.h"

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
   if(fabs(s-1.0)>=0.0001) printf("%10s has allele frequencies summed up to %f \n",m.name,s);   
   for(i=0;i<m.n;++i)
   {
      m.freq[i]/=s;
   }
   locus[j]=m;++j;
}
fclose(fp);
nloci=j; 

return 0;

}

int main(int argc, char *argv[])
{
char frqname[20],freqname[20],mrkname[20];
FILE *header,*name;
int i,j,k;

for(k=1;k<=23;k++)
{
   if(k<10) {sprintf(frqname,"ms0%d.frq",k);sprintf(freqname,"ms0%d.freq",k);sprintf(mrkname,"header0%d",k);}
   else {sprintf(frqname,"ms%d.frq",k);sprintf(freqname,"ms%d.freq",k);sprintf(mrkname,"header%d",k);}
   getfrq(frqname);
   header=fopen(mrkname,"w");
   name=fopen(freqname,"w");
   if(!name||!header)
   {
     fprintf(stderr,"error to open file for chromosome %d\n",k);
     return 1;
   }
   fprintf(header,"id,");
   for(i=0;i<nloci-1;i++) fprintf(header,"%s,",locus[i].name);
   fprintf(header,"%s\n",locus[nloci-1].name);
   fclose(header);
   for(i=0;i<nloci;i++) 
   {
     fprintf(name,"%s ",locus[i].name);
     for(j=0;j<locus[i].n;j++) fprintf(name,"%d %f ",locus[i].size[j],locus[i].freq[j]);
     fprintf(name,"\n");
   }
   fclose(name);
}
return 0;
}
