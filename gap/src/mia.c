#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cline.h"
#include "hap.h"

#define version 1.0
#define MAX_NAME_LEN 30

typedef struct {
  int mimp,chr;
  char id[MAX_NAME_LEN],**namei;
  double prob;
} so_def;
so_def *so_list,*so_list_t;
static int n_loci;

void qsorts(long int,long int);

void mia_c(char **hapfile,char **assfile, char **miafile, int *so, int *ns, int *mi, int *allsnps, int *sas)
{
  FILE *infile,*outfile,*sasfile;
  char ofname[MAX_FILENAME_LEN];
  char tpname[MAX_FILENAME_LEN];
  char line[MAX_LINE_LEN], rol[MAX_LINE_LEN];
  char namei[MAX_NAME_LEN],a[3],**names,***hapnames;
  int i,j,k,l,n_hap,mimp;
  double **haptable;
  double freq;
  double hap_mean,hap_sd,hap_min,hap_max;
  char id[MAX_NAME_LEN];

  mimp=*mi;
  if(mimp==0) {
    REprintf("\nPlease specify -mi #, # being number of imputations\n");
    return;
  }

  sprintf(ofname,"%s",*miafile);
  if(!ofname[0]) strcpy(ofname,"mia.out");
  outfile=fopen(ofname,"w");
  if(!outfile) goto open_error;

  n_loci=0;
  n_hap=0;

  sasfile=NULL;
  if(!*so) /*haplotype order*/
  {
    sprintf(tpname,"%s.%03d",*hapfile,1);
    infile=fopen(tpname,"r");
    if(!infile) goto read_error;
    else {
      fgets(line, MAX_LINE_LEN, infile);
      fprintf(outfile,"Haplotypes\n\n   ID");
      while(sscanf(line," %s %[^\n]",namei,rol)>1) {
        strcpy(line,rol);
        fprintf(outfile," %s",namei);
        n_loci++;
      }
      fprintf(outfile,"\n\n");
      while(fgets(line,MAX_LINE_LEN,infile)) n_hap++;
      rewind(infile);
    }
    hapnames=(char***)malloc(n_hap*sizeof(char**));
    if(!hapnames) goto no_room;
    for(i=0;i<n_hap;i++) {
      hapnames[i]=(char**)malloc(n_loci*sizeof(char*));
      if(!hapnames[i]) goto no_room;
      for(j=0;j<n_loci;j++) {
        hapnames[i][j]=(char*)malloc(3*sizeof(char));
        if(!hapnames[i][j]) goto no_room;
      }
      for(j=0;j<n_loci;j++) strcpy(hapnames[i][j]," ");
    }
    names=(char**)malloc(n_loci*sizeof(char*));
    if(!names) goto no_room;
    for(i=0;i<n_loci;i++) {
      names[i]=(char*)malloc(MAX_NAME_LEN*sizeof(char));
      if(!names[i]) goto no_room;
    }
    haptable=(double**) malloc(n_hap*sizeof(double*));
    if(!haptable) goto no_room;
    for(i=0;i<n_hap;i++) {
      haptable[i]=(double *)malloc(mimp*sizeof(double));
      if(!haptable[i]) goto no_room;
      for(j=0;j<mimp;j++) haptable[i][j]=0;
    }
    for(i=1;i<=mimp;i++) {
      if(i>1) {
        sprintf(tpname,"%s.%03d",*hapfile,i);
        infile=fopen(tpname,"r");
      }
      if(!infile) goto read_error;
      else {
        fgets(line, MAX_LINE_LEN, infile);
        j=0;
        if(i==1)
          while(sscanf(line," %s %[^\n]",namei,rol)>1) {
            strcpy(names[j],namei);
            strcpy(line,rol);
            ++j;
          }
        k=0;
        while(fgets(line, MAX_LINE_LEN, infile)) {
          if(i==1) fprintf(outfile,"%5d ",k+1);
          for(j=0;j<n_loci;j++) {
            sscanf(line,"%s %[^\n]",a,rol);
            if(i==1) {
              strcpy(hapnames[k][j],a);
              if(*allsnps) fprintf(outfile,"%1s",a);
              else fprintf(outfile,"%2s",a);
            }
            strcpy(line,rol);
          }
          if(i==1) fprintf(outfile,"\n");
          sscanf(rol,"%lf",&freq);
          haptable[k][i-1]=freq;
          k++;
        }
        fclose(infile);
      }
    }
    fprintf(outfile,"\nSummary statistics of haplotype frequencies\n");
    fprintf(outfile,"\n   ID     Mean     s.d.      Min      Max\n\n");
    for(i=0;i<n_hap;i++) {
      hap_mean=0;
      hap_sd=0;
      hap_min=1.0;
      hap_max=0.0;
      for(j=0;j<mimp;j++) {
        freq=haptable[i][j];
        hap_mean+=freq;
        hap_sd+=freq*freq;
        if(freq<hap_min) hap_min=freq;
        if(freq>hap_max) hap_max=freq;
      }
      if(mimp>1) hap_sd=sqrt((hap_sd-hap_mean*hap_mean/mimp)/(mimp-1));
      hap_mean/=mimp;
      fprintf(outfile,"%5d %f %f",i+1,hap_mean,hap_sd);
      fprintf(outfile," %f %f\n",hap_min,hap_max);
    }
    fprintf(outfile,"\n\nThe actual haplotype frequencies by imputations\n\n");
    fprintf(outfile,"   ID");
    for(j=0;j<mimp;j++) fprintf(outfile," %8d",j+1);
    fprintf(outfile,"\n\n");
    if(*sas) {
      sprintf(tpname,"%s.sas",*hapfile);
      sasfile=fopen(tpname,"w");
      if(!sasfile) goto open_error;
      fprintf(sasfile,"/*produced from HAP output by JH Zhao*/\n\n");
      fprintf(sasfile,"data hapfreq;\n");
      fprintf(sasfile,"input id");
      k=0;
      for(i=0;i<n_loci;i++) {
        if(strcmp(names[i],"0")>=0&&strcmp(names[i],"9")<=0)
          fprintf(sasfile," site%-d$",++k);
        else fprintf(sasfile," %s$",names[i]);
      }
      if(mimp==1) fprintf(sasfile," h1");
      else fprintf(sasfile," h1-h%-d",mimp);
      fprintf(sasfile,";\n");
      fprintf(sasfile,"cards;\n");
    }
    for(i=0;i<n_hap;i++) {
      fprintf(outfile,"%5d",i+1);
      for(j=0;j<mimp;j++) fprintf(outfile," %f",haptable[i][j]);
      fprintf(outfile,"\n");
      if(*sas) {
        fprintf(sasfile,"%5d",i+1);
        for(j=0;j<n_loci;j++) fprintf(sasfile,"%3s",hapnames[i][j]);
        for(j=0;j<mimp;j++) fprintf(sasfile," %f",haptable[i][j]);
        fprintf(sasfile,"\n");
      }
    }
    if(*sas) {
      fprintf(sasfile,";\n");
      REprintf("\nSAS program file has been written to %s\n",tpname);
      fclose(sasfile);
    }
    for(i=0;i<n_hap;i++) {
      for(j=0;j<n_loci;j++) free(hapnames[i][j]);
      free(hapnames[i]);
    }
    free(hapnames);
    for(i=0;i<n_hap;i++) free(haptable[i]);
    free(haptable);
  }
  else  /*subject order*/
  {
    fprintf(outfile,"Site names\n\n");
    for(i=1;i<=mimp;i++) {
      sprintf(tpname,"%s.%03d",*assfile,i);
      infile=fopen(tpname,"r");
      if(!infile) goto read_error;
      else {
        fgets(line, MAX_LINE_LEN, infile);
        if(i==1) {
          sscanf(line, "%*s %*s %[^\n]",rol);
          strcpy(line,rol);
          while(sscanf(line," %s %[^\n]",namei,rol)>1) {
            strcpy(line,rol);
            fprintf(outfile," %s",namei);
            n_loci++;
          }
          fprintf(outfile,"\n");
        }
        while(fgets(line,MAX_LINE_LEN,infile)) n_hap++;
        fclose(infile);
      }
    }
    names=(char**)malloc(n_loci*sizeof(char*));
    if(!names) goto no_room;
    for(i=0;i<n_loci;i++) {
      names[i]=(char*)malloc(MAX_NAME_LEN*sizeof(char));
      if(!names[i]) goto no_room;
    }
    so_list=(so_def*)malloc(n_hap*sizeof(so_def));
    if(!so_list) goto no_room;
    for(i=0;i<n_hap;i++) {
      so_list[i].namei=(char**)malloc(n_loci*sizeof(char*));
      if(!so_list[i].namei) goto no_room;
      for(j=0;j<n_loci;j++) {
        so_list[i].namei[j]=(char*)malloc(3*sizeof(char));
        if(!so_list[i].namei[j]) goto no_room;
      }
      for(j=0;j<n_loci;j++) strcpy(so_list[i].namei[j]," ");
    }
    so_list_t=(so_def*)malloc(sizeof(so_def));
    if(!so_list_t) goto no_room;
    so_list_t[0].namei=(char**)malloc(n_loci*sizeof(char*));
    if(!so_list_t[0].namei) goto no_room;
    for(j=0;j<n_loci;j++) {
      so_list_t[0].namei[j]=(char*)malloc(3*sizeof(char));
      if(!so_list_t[0].namei[j]) goto no_room;
    }
    l=0;
    for(i=1;i<=mimp;i++) {
      sprintf(tpname,"%s.%03d",*assfile,i);
      infile=fopen(tpname,"r");
      if(infile) {
        fgets(line,MAX_LINE_LEN,infile);
        j=0;
        if(i==1) sscanf(line, "%*s %*s %[^\n]",rol);
        strcpy(line,rol);
        while(sscanf(line," %s %[^\n]",namei,rol)>1) {
            strcpy(names[j],namei);
            strcpy(line,rol);
            ++j;
        }
        while(fgets(line, MAX_LINE_LEN, infile)&&
              sscanf(line,"%s %d %[^\n]",id,&k,rol)) {
          so_list[l].mimp=i;
          strcpy(so_list[l].id,id);
          so_list[l].chr=k;
          strcpy(line,rol);
          for(j=0;j<n_loci;j++) {
            sscanf(line,"%s %[^\n]",a,rol);
            strcpy(so_list[l].namei[j],a);
            strcpy(line,rol);
          }
          sscanf(rol,"%lf",&so_list[l].prob);
          ++l;
        }
        fclose(infile);
      } else goto read_error;
    }
    if(*sas) {
      sprintf(tpname,"%s.sas",*assfile);
      sasfile=fopen(tpname,"w");
      if(!sasfile) goto open_error;
      fprintf(sasfile,"/*produced from HAP output by JH Zhao*/\n\n");
      fprintf(sasfile,"data subject;\n");
      fprintf(sasfile,"input id$ imp prob chr");
      k=0;
      for(j=0;j<n_loci;j++) {
        if(strcmp(names[j],"0")>=0&&strcmp(names[j],"9")<=0)
          fprintf(sasfile," site%-d$",++k);
        else fprintf(sasfile," %s$",names[j]);
      }
      fprintf(sasfile,";\n");
      fprintf(sasfile,"cards;\n");
      for(i=0;i<n_hap;i++) {
        fprintf(sasfile," %8s",so_list[i].id);
        fprintf(sasfile," %3d",so_list[i].mimp);
        fprintf(sasfile," %f",so_list[i].prob);
        fprintf(sasfile," %1d ",so_list[i].chr);
        for(j=0;j<n_loci;j++) fprintf(sasfile," %2s",so_list[i].namei[j]);
        fprintf(sasfile,"\n");
      }
      fprintf(sasfile,";\n");
      REprintf("\nSAS program has been written to %s\n",tpname);
      fclose(sasfile);
    }
    if(!*ns) qsorts(0,n_hap-1);
    fprintf(outfile,"\nAssignment of haplotypes by ID & imputation\n");
    fprintf(outfile,"\n       ID imp  prob  chr haplotype\n\n");
    for(i=0;i<n_hap;i++) {
      fprintf(outfile," %8s",so_list[i].id);
      fprintf(outfile," %3d",so_list[i].mimp);
      fprintf(outfile," %f",so_list[i].prob);
      fprintf(outfile," %1d ",so_list[i].chr);
      for(j=0;j<n_loci;j++) {
        if(*allsnps) fprintf(outfile,"%1s",so_list[i].namei[j]);
        else fprintf(outfile,"%2s",so_list[i].namei[j]);
      }
      fprintf(outfile,"\n");
    }
    for(i=0;i<n_hap;i++) {
      for(j=0;j<n_loci;j++) free(so_list[i].namei[j]);
      free(so_list[i].namei);
    }
    free(so_list);
    free(so_list_t);
  }
  for(j=0;j<n_loci;j++) free(names[j]);
  free(names);
  fclose(outfile);
  REprintf("Output has been written to %s",ofname);

  return;

  open_error:
  REprintf("Error opening file\n");
  return;

  read_error:
  REprintf("Error reading file %s\n",tpname);
  return;

  no_room:
  REprintf("Insufficient memory\n");
  return;

}

#ifdef executable
int main (int argc, char **argv)
{
  FILE *infile,*outfile,*sasfile;
  char ifname[MAX_FILENAME_LEN], ofname[MAX_FILENAME_LEN];
  char tpname[MAX_FILENAME_LEN];
  char line[MAX_LINE_LEN], rol[MAX_LINE_LEN];
  char namei[MAX_NAME_LEN],a[3],**names,***hapnames;
  int i,j,k,l,n_hap;
  int so=0,ns=0,mimp=0,allsnps=0,sas=0;
  double **haptable;
  double freq;
  double hap_mean,hap_sd,hap_min,hap_max;
  char id[MAX_NAME_LEN];

  if (
      get_flag(argc, argv, "so",0 , &so) < 0 ||
      get_flag(argc, argv, "ns",0 , &ns) < 0 ||
      get_flag(argc, argv, "mi",1 , &mimp) < 0 ||
      get_flag(argc, argv, "as",0 , &allsnps) < 0 ||
      get_flag(argc, argv, "sas",0 , &sas) < 0 ||
      !get_arg(argc, argv, ifname)
     ) {
    REprintf("\nMIANALYZE version %.2f JH Zhao 2002\n\n",version);
    REprintf("Usage: %s [options]", argv[0]);
    REprintf(" input-file-root [output-file]\n");
    REprintf("\n Where options (defaults) are:\n\n");
    REprintf("\t-so \tTally haplotypes by subject order\n");
    REprintf("\t-ns \tDo not sort by individual ID\n");
    REprintf("\t-mi #\tNumber of imputations use in HAP\n");
    REprintf("\t-as \tAll markers are SNPs\n");
    REprintf("\t-sas \tTo generate SAS data step statements\n");
    REprintf("\n input-file-root is output filename");
    REprintf(" produced by HAP with -mi# and -ss options\n");
    Rprintf("\n");
    Rprintf("Max. length of finenames   = %d\n",MAX_FILENAME_LEN);
    Rprintf("Max. length of input lines = %d\n",MAX_LINE_LEN);
    Rprintf("Max. length of locus names = %d\n",MAX_NAME_LEN);
    Rprintf("\n");
    return 1;
  }

  if(mimp==0) {
    REprintf("\nPlease specify -mi #, # being number of imputations\n");
    return 1;
  }
  get_arg(argc, argv, ofname);

  if(!ofname[0]) strcpy(ofname,"mia.out");

  outfile=fopen(ofname,"w");
  if(!outfile) goto open_error;

  n_loci=0;
  n_hap=0;

  if(!so) /*haplotype order*/
  {
    sprintf(tpname,"%s.%03d",ifname,1);
    infile=fopen(tpname,"r");
    if(!infile) goto read_error;
    else {
      fgets(line, MAX_LINE_LEN, infile);
      fprintf(outfile,"Haplotypes\n\n   ID");
      while(sscanf(line," %s %[^\n]",namei,rol)>1) {
        strcpy(line,rol);
        fprintf(outfile," %s",namei);
        n_loci++;
      }
      fprintf(outfile,"\n\n");
      while(fgets(line,MAX_LINE_LEN,infile)) n_hap++;
      rewind(infile);
    }
    hapnames=(char***)malloc(n_hap*sizeof(char**));
    if(!hapnames) goto no_room;
    for(i=0;i<n_hap;i++) {
      hapnames[i]=(char**)malloc(n_loci*sizeof(char*));
      if(!hapnames[i]) goto no_room;
      for(j=0;j<n_loci;j++) {
        hapnames[i][j]=(char*)malloc(3*sizeof(char));
        if(!hapnames[i][j]) goto no_room;
      }
      for(j=0;j<n_loci;j++) strcpy(hapnames[i][j]," ");
    }
    names=(char**)malloc(n_loci*sizeof(char*));
    if(!names) goto no_room;
    for(i=0;i<n_loci;i++) {
      names[i]=(char*)malloc(MAX_NAME_LEN*sizeof(char));
      if(!names[i]) goto no_room;
    }
    haptable=(double**) malloc(n_hap*sizeof(double*));
    if(!haptable) goto no_room;
    for(i=0;i<n_hap;i++) {
      haptable[i]=(double *)malloc(mimp*sizeof(double));
      if(!haptable[i]) goto no_room;
      for(j=0;j<mimp;j++) haptable[i][j]=0;
    }
    for(i=1;i<=mimp;i++) {
      if(i>1) {
        sprintf(tpname,"%s.%03d",ifname,i);
        infile=fopen(tpname,"r");
      }
      if(!infile) goto read_error;
      else {
        fgets(line, MAX_LINE_LEN, infile);
        j=0;
        if(i==1)
          while(sscanf(line," %s %[^\n]",namei,rol)>1) {
            strcpy(names[j],namei);
            strcpy(line,rol);
            ++j;
          }
        k=0;
        while(fgets(line, MAX_LINE_LEN, infile)) {
          if(i==1) fprintf(outfile,"%5d ",k+1);
          for(j=0;j<n_loci;j++) {
            sscanf(line,"%s %[^\n]",a,rol);
            if(i==1) {
              strcpy(hapnames[k][j],a);
              if(allsnps) fprintf(outfile,"%1s",a);
              else fprintf(outfile,"%2s",a);
            }
            strcpy(line,rol);
          }
          if(i==1) fprintf(outfile,"\n");
          sscanf(rol,"%f",&freq);
          haptable[k][i-1]=freq;
          k++;
        }
        fclose(infile);
      }
    }
    fprintf(outfile,"\nSummary statistics of haplotype frequencies\n");
    fprintf(outfile,"\n   ID     Mean     s.d.      Min      Max\n\n");
    for(i=0;i<n_hap;i++) {
      hap_mean=0;
      hap_sd=0;
      hap_min=1.0;
      hap_max=0.0;
      for(j=0;j<mimp;j++) {
        freq=haptable[i][j];
        hap_mean+=freq;
        hap_sd+=freq*freq;
        if(freq<hap_min) hap_min=freq;
        if(freq>hap_max) hap_max=freq;
      }
      if(mimp>1) hap_sd=sqrt((hap_sd-hap_mean*hap_mean/mimp)/(mimp-1));
      hap_mean/=mimp;
      fprintf(outfile,"%5d %f %f",i+1,hap_mean,hap_sd);
      fprintf(outfile," %f %f\n",hap_min,hap_max);
    }
    fprintf(outfile,"\n\nThe actual haplotype frequencies by imputations\n\n");
    fprintf(outfile,"   ID");
    for(j=0;j<mimp;j++) fprintf(outfile," %8d",j+1);
    fprintf(outfile,"\n\n");
    if(sas) {
      sprintf(tpname,"%s.sas",ifname);
      sasfile=fopen(tpname,"w");
      if(!sasfile) goto open_error;
      fprintf(sasfile,"/*produced from HAP output by JH Zhao*/\n\n");
      fprintf(sasfile,"data hapfreq;\n");
      fprintf(sasfile,"input id");
      k=0;
      for(i=0;i<n_loci;i++) {
        if(strcmp(names[i],"0")>=0&&strcmp(names[i],"9")<=0)
          fprintf(sasfile," site%-d$",++k);
        else fprintf(sasfile," %s$",names[i]);
      }
      if(mimp==1) fprintf(sasfile," h1");
      else fprintf(sasfile," h1-h%-d",mimp);
      fprintf(sasfile,";\n");
      fprintf(sasfile,"cards;\n");
    }
    for(i=0;i<n_hap;i++) {
      fprintf(outfile,"%5d",i+1);
      for(j=0;j<mimp;j++) fprintf(outfile," %f",haptable[i][j]);
      fprintf(outfile,"\n");
      if(sas) {
        fprintf(sasfile,"%5d",i+1);
        for(j=0;j<n_loci;j++) fprintf(sasfile,"%3s",hapnames[i][j]);
        for(j=0;j<mimp;j++) fprintf(sasfile," %f",haptable[i][j]);
        fprintf(sasfile,"\n");
      }
    }
    if(sas) {
      fprintf(sasfile,";\n");
      REprintf("\nSAS program file has been written to %s\n",tpname);
      fclose(sasfile);
    }
    for(i=0;i<n_hap;i++) {
      for(j=0;j<n_loci;j++) free(hapnames[i][j]);
      free(hapnames[i]);
    }
    free(hapnames);
    for(i=0;i<n_hap;i++) free(haptable[i]);
    free(haptable);
  }
  else  /*subject order*/
  {
    fprintf(outfile,"Site names\n\n");
    for(i=1;i<=mimp;i++) {
      sprintf(tpname,"%s.%03d",ifname,i);
      infile=fopen(tpname,"r");
      if(!infile) goto read_error;
      else {
        fgets(line, MAX_LINE_LEN, infile);
        if(i==1) {
          sscanf(line, "%*s %*s %[^\n]",rol);
          strcpy(line,rol);
          while(sscanf(line," %s %[^\n]",namei,rol)>1) {
            strcpy(line,rol);
            fprintf(outfile," %s",namei);
            n_loci++;
          }
          fprintf(outfile,"\n");
        }
        while(fgets(line,MAX_LINE_LEN,infile)) n_hap++;
        fclose(infile);
      }
    }
    names=(char**)malloc(n_loci*sizeof(char*));
    if(!names) goto no_room;
    for(i=0;i<n_loci;i++) {
      names[i]=(char*)malloc(MAX_NAME_LEN*sizeof(char));
      if(!names[i]) goto no_room;
    }
    so_list=(so_def*)malloc(n_hap*sizeof(so_def));
    if(!so_list) goto no_room;
    for(i=0;i<n_hap;i++) {
      so_list[i].namei=(char**)malloc(n_loci*sizeof(char*));
      if(!so_list[i].namei) goto no_room;
      for(j=0;j<n_loci;j++) {
        so_list[i].namei[j]=(char*)malloc(3*sizeof(char));
        if(!so_list[i].namei[j]) goto no_room;
      }
      for(j=0;j<n_loci;j++) strcpy(so_list[i].namei[j]," ");
    }
    so_list_t=(so_def*)malloc(sizeof(so_def));
    if(!so_list_t) goto no_room;
    so_list_t[0].namei=(char**)malloc(n_loci*sizeof(char*));
    if(!so_list_t[0].namei) goto no_room;
    for(j=0;j<n_loci;j++) {
      so_list_t[0].namei[j]=(char*)malloc(3*sizeof(char));
      if(!so_list_t[0].namei[j]) goto no_room;
    }
    l=0;
    for(i=1;i<=mimp;i++) {
      sprintf(tpname,"%s.%03d",ifname,i);
      infile=fopen(tpname,"r");
      if(infile) {
        fgets(line,MAX_LINE_LEN,infile);
        j=0;
        if(i==1)
          sscanf(line, "%*s %*s %[^\n]",rol);
          strcpy(line,rol);
          while(sscanf(line," %s %[^\n]",namei,rol)>1) {
            strcpy(names[j],namei);
            strcpy(line,rol);
            ++j;
          }
        while(fgets(line, MAX_LINE_LEN, infile)&&
              sscanf(line,"%s %d %[^\n]",id,&k,rol)) {
          so_list[l].mimp=i;
          strcpy(so_list[l].id,id);
          so_list[l].chr=k;
          strcpy(line,rol);
          for(j=0;j<n_loci;j++) {
            sscanf(line,"%s %[^\n]",a,rol);
            strcpy(so_list[l].namei[j],a);
            strcpy(line,rol);
          }
          sscanf(rol,"%f",&so_list[l].prob);
          ++l;
        }
        fclose(infile);
      } else goto read_error;
    }
    if(sas) {
      sprintf(tpname,"%s.sas",ifname);
      sasfile=fopen(tpname,"w");
      if(!sasfile) goto open_error;
      fprintf(sasfile,"/*produced from HAP output by JH Zhao*/\n\n");
      fprintf(sasfile,"data subject;\n");
      fprintf(sasfile,"input id$ imp prob chr");
      k=0;
      for(j=0;j<n_loci;j++) {
        if(strcmp(names[j],"0")>=0&&strcmp(names[j],"9")<=0)
          fprintf(sasfile," site%-d$",++k);
        else fprintf(sasfile," %s$",names[j]);
      }
      fprintf(sasfile,";\n");
      fprintf(sasfile,"cards;\n");
      for(i=0;i<n_hap;i++) {
        fprintf(sasfile," %8s",so_list[i].id);
        fprintf(sasfile," %3d",so_list[i].mimp);
        fprintf(sasfile," %f",so_list[i].prob);
        fprintf(sasfile," %1d ",so_list[i].chr);
        for(j=0;j<n_loci;j++) fprintf(sasfile," %2s",so_list[i].namei[j]);
        fprintf(sasfile,"\n");
      }
      fprintf(sasfile,";\n");
      REprintf("\nSAS program has been written to %s\n",tpname);
      fclose(sasfile);
    }
    if(!ns) qsorts(0,n_hap-1);
    fprintf(outfile,"\nAssignment of haplotypes by ID & imputation\n");
    fprintf(outfile,"\n       ID imp  prob  chr haplotype\n\n");
    for(i=0;i<n_hap;i++) {
      fprintf(outfile," %8s",so_list[i].id);
      fprintf(outfile," %3d",so_list[i].mimp);
      fprintf(outfile," %f",so_list[i].prob);
      fprintf(outfile," %1d ",so_list[i].chr);
      for(j=0;j<n_loci;j++) {
        if(allsnps) fprintf(outfile,"%1s",so_list[i].namei[j]);
        else fprintf(outfile,"%2s",so_list[i].namei[j]);
      }
      fprintf(outfile,"\n");
    }
    for(i=0;i<n_hap;i++) {
      for(j=0;j<n_loci;j++) free(so_list[i].namei[j]);
      free(so_list[i].namei);
    }
    free(so_list);
    free(so_list_t);
  }
  for(j=0;j<n_loci;j++) free(names[j]);
  free(names);
  fclose(outfile);
  REprintf("Output has been written to %s",ofname);

  return 0;

  open_error:
  REprintf("Error opening file\n");
  return 1;

  read_error:
  REprintf("Error reading file\n");
  return 1;

  no_room:
  REprintf("Insufficient memory\n");
  return 1;

}
#endif

void qsorts(long int from,long int to)
/************************************************/
/* David Brunskill and John Turner (1996)       */
/* Understanding Algorithms and Data Structures */
/* McGraw-Hill                                  */
/************************************************/
{
long int i,pivot;

if(to>from)
{
  pivot=from;
  for(i=from+1;i<=to;++i)
  {
    so_list_t[0]=so_list[i];
    if(strcmp(so_list_t[0].id,so_list[pivot].id)<=0)
    {
      so_list[i]=so_list[pivot+1];
      so_list[pivot+1]=so_list[pivot];
      so_list[pivot]=so_list_t[0];
      pivot++;
    }
  }
  qsorts(from,pivot-1);
  qsorts(pivot+1,to);
}
}
